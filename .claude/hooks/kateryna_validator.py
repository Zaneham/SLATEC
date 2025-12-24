#!/usr/bin/env python3
"""
Kateryna epistemic validator hook for Claude Code.

Validates Claude's outputs against RAG sources to detect potential
hallucinations (confident responses with weak grounding).

Epistemic scoring:
  +1: Confident + Strong grounding = Trust
   0: Uncertain + Weak grounding   = Appropriate uncertainty
  -1: Confident + Weak grounding   = HALLUCINATION RISK

Now with Hopper knowledge base integration for real RAG validation.
"""
import json
import re
import sqlite3
import sys
import os
from pathlib import Path
from typing import Tuple, List, Optional

try:
    from kateryna import EpistemicDetector, TernaryState
    KATERYNA_AVAILABLE = True
except ImportError:
    KATERYNA_AVAILABLE = False

# Hopper knowledge base path
HOPPER_DB_PATH = Path(os.environ.get("HOPPER_DB_PATH", "C:/dev/hopper/hopper_index.db"))

# Stop words for term extraction
STOP_WORDS = {
    'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
    'of', 'with', 'by', 'from', 'is', 'are', 'was', 'were', 'be', 'been',
    'being', 'have', 'has', 'had', 'do', 'does', 'did', 'will', 'would',
    'could', 'should', 'may', 'might', 'must', 'shall', 'can', 'this',
    'that', 'these', 'those', 'it', 'its', 'if', 'then', 'else', 'when',
    'where', 'which', 'who', 'what', 'how', 'why', 'all', 'each', 'every',
    'both', 'few', 'more', 'most', 'other', 'some', 'such', 'no', 'not',
    'only', 'own', 'same', 'so', 'than', 'too', 'very', 'just', 'also',
    'now', 'here', 'there', 'use', 'used', 'using', 'end', 'return',
    'function', 'subroutine', 'call', 'intent', 'real', 'integer', 'type',
}


class HopperSearch:
    """Query Hopper knowledge base for RAG validation."""

    def __init__(self, db_path: Path = HOPPER_DB_PATH):
        self.db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None

    def _get_conn(self) -> sqlite3.Connection:
        """Lazy connection to database."""
        if self._conn is None:
            if not self.db_path.exists():
                raise FileNotFoundError(f"Hopper database not found: {self.db_path}")
            self._conn = sqlite3.connect(str(self.db_path))
        return self._conn

    def extract_terms(self, content: str, max_terms: int = 10) -> List[Tuple[str, int]]:
        """
        Extract key search terms from content with importance scores.

        Returns list of (term, importance) tuples, sorted by importance.
        Higher importance = term looks more "technical" or "specific".
        """
        # Find potential technical terms (capitalised, underscored, or longer words)
        words = re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', content)

        # Filter and score terms
        term_scores = {}
        for word in words:
            lower = word.lower()
            if lower in STOP_WORDS or len(word) < 3:
                continue

            # Prioritise:
            # - ALL CAPS (likely acronyms/constants): +3
            # - CamelCase or with underscores: +2
            # - Longer words: +1
            score = 1  # Base score
            if word.isupper() and len(word) >= 2:
                score += 3
            if '_' in word or (word[0].isupper() and not word.isupper()):
                score += 2
            if len(word) >= 6:
                score += 1

            term_scores[lower] = max(term_scores.get(lower, 0), score)

        # Sort by score descending, take top N
        sorted_terms = sorted(term_scores.items(), key=lambda x: -x[1])
        return sorted_terms[:max_terms]

    def search(self, terms: List[Tuple[str, int]], limit: int = 20) -> Tuple[float, int]:
        """
        Search Hopper FTS5 index for terms with importance weighting.

        Key insight: unfound *important* terms should REDUCE confidence.
        A fabricated routine name like "XYZFOOBAR" (high importance, all caps)
        not being in the KB is a red flag, even if common terms match.

        Returns:
            (retrieval_confidence, chunks_found)
            - retrieval_confidence: 0.0-1.0, penalised for unfound important terms
            - chunks_found: total number of matching chunks
        """
        if not terms:
            return 0.0, 0

        conn = self._get_conn()
        cursor = conn.cursor()

        total_chunks = 0
        weighted_found = 0.0
        weighted_unfound = 0.0
        total_weight = 0.0

        for term, importance in terms:
            total_weight += importance
            try:
                # Count total matches for this term (to detect common words)
                cursor.execute("""
                    SELECT COUNT(*) FROM document_fts WHERE document_fts MATCH ?
                """, (term,))
                doc_count = cursor.fetchone()[0]

                cursor.execute("""
                    SELECT COUNT(*) FROM snippet_fts WHERE snippet_fts MATCH ?
                """, (term,))
                snip_count = cursor.fetchone()[0]

                term_matches = doc_count + snip_count

                # If term matches too many docs, it's probably a common word
                # that doesn't provide real grounding (e.g. "correct", "definitely")
                # Max 500 matches = specific enough to count as grounding
                if term_matches > 0 and term_matches <= 500:
                    weighted_found += importance
                    total_chunks += term_matches
                elif term_matches == 0:
                    # Unfound term - penalise by importance
                    weighted_unfound += importance
                # else: too common, don't count either way

            except sqlite3.OperationalError:
                # Skip malformed terms
                continue

        if total_weight == 0:
            return 0.0, 0

        # Confidence = (found weight - unfound weight) / total weight
        # This means unfound important terms actively reduce confidence
        # Range: -1.0 to 1.0, then normalised to 0.0-1.0
        raw_score = (weighted_found - weighted_unfound) / total_weight
        retrieval_confidence = (raw_score + 1.0) / 2.0  # Normalise to 0-1

        return retrieval_confidence, total_chunks

    def close(self):
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None


def validate_content(content: str, context: str = "") -> dict:
    """
    Validate content using Kateryna's epistemic detector with Hopper RAG.

    Args:
        content: The text Claude produced
        context: Optional context about what was being asked

    Returns:
        dict with 'decision' and optional 'reason'
    """
    if not KATERYNA_AVAILABLE:
        # Kateryna not installed - allow but warn
        return {
            "decision": "approve",
            "systemMessage": "Kateryna not installed - skipping epistemic validation"
        }

    # Query Hopper for RAG context
    retrieval_confidence = float(os.environ.get("KATERYNA_RETRIEVAL_CONFIDENCE", "0.5"))
    chunks_found = int(os.environ.get("KATERYNA_CHUNKS_FOUND", "1"))

    hopper = None
    try:
        if HOPPER_DB_PATH.exists():
            hopper = HopperSearch()
            terms = hopper.extract_terms(content)
            if terms:
                retrieval_confidence, chunks_found = hopper.search(terms)
    except Exception as e:
        # Hopper unavailable - fall back to env vars
        pass
    finally:
        if hopper:
            hopper.close()

    detector = EpistemicDetector()

    # Analyse the content with real RAG values
    state = detector.analyze(
        text=content,
        question=context or "Validate this response",
        retrieval_confidence=retrieval_confidence,
        chunks_found=chunks_found
    )

    if state.is_danger_zone:
        # -1: Confident bullshit detected
        return {
            "decision": "block",
            "reason": f"Kateryna epistemic check failed: {state.reason}",
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": (
                    f"Epistemic state: {state.state.name} - "
                    f"RAG confidence: {retrieval_confidence:.2f}, chunks: {chunks_found}. "
                    "Consider adding uncertainty or verifying against sources."
                )
            }
        }

    # +1 or 0: Acceptable
    return {"decision": "approve"}


def main():
    """Main hook handler - receives JSON via stdin."""
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError:
        # No valid input - allow
        sys.exit(0)

    # Only validate PostToolUse events
    if input_data.get("hook_event_name") != "PostToolUse":
        sys.exit(0)

    tool_name = input_data.get("tool_name", "")
    tool_input = input_data.get("tool_input", {})

    # Extract content based on tool type
    content = None
    if tool_name == "Write":
        content = tool_input.get("content", "")
    elif tool_name == "Edit":
        content = tool_input.get("new_string", "")

    # Skip if no content to validate
    if not content or len(content) < 50:  # Skip trivial edits
        print(json.dumps({"decision": "approve"}))
        sys.exit(0)

    # Run validation
    result = validate_content(content)
    print(json.dumps(result))
    sys.exit(0)


if __name__ == "__main__":
    main()
