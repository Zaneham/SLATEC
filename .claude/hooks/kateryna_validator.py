#!/usr/bin/env python3
"""
Kateryna epistemic validator hook for Claude Code.

Validates Claude's outputs against RAG sources to detect potential
hallucinations (confident responses with weak grounding).

Epistemic scoring:
  +1: Confident + Strong grounding = Trust
   0: Uncertain + Weak grounding   = Appropriate uncertainty
  -1: Confident + Weak grounding   = HALLUCINATION RISK
"""
import json
import sys
import os

try:
    from kateryna import EpistemicDetector, TernaryState
    KATERYNA_AVAILABLE = True
except ImportError:
    KATERYNA_AVAILABLE = False


def validate_content(content: str, context: str = "") -> dict:
    """
    Validate content using Kateryna's epistemic detector.

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

    detector = EpistemicDetector()

    # Analyse the content
    # TODO: Connect to your actual RAG system to get retrieval_confidence and chunks_found
    # For now, using defaults that assume no RAG context available
    state = detector.analyze(
        text=content,
        question=context or "Validate this response",
        retrieval_confidence=float(os.environ.get("KATERYNA_RETRIEVAL_CONFIDENCE", "0.5")),
        chunks_found=int(os.environ.get("KATERYNA_CHUNKS_FOUND", "1"))
    )

    if state.is_danger_zone:
        # -1: Confident bullshit detected
        return {
            "decision": "block",
            "reason": f"Kateryna epistemic check failed: {state.reason}",
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": f"Epistemic state: {state.value} - Consider adding uncertainty or verifying against sources."
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
