#!/usr/bin/env python
"""Find files with active GOTO statements (not in comments)."""

import os
import re

root = 'C:/dev/slatec-modern/src/modern'
goto_pattern = re.compile(r'^\s*(IF.*)?GOTO\s+\d+', re.IGNORECASE)

files_by_count = {}

for dirpath, dirnames, filenames in os.walk(root):
    for filename in filenames:
        if filename.endswith('.f90'):
            filepath = os.path.join(dirpath, filename)
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            count = 0
            for line in lines:
                # Skip comment lines
                stripped = line.strip()
                if stripped.startswith('!'):
                    continue
                # Check if line has GOTO (not in a comment part)
                if '!' in line:
                    code_part = line.split('!')[0]
                else:
                    code_part = line
                if goto_pattern.search(code_part):
                    count += 1

            if count > 0:
                if count not in files_by_count:
                    files_by_count[count] = []
                files_by_count[count].append(filepath.replace(root + '/', ''))

# Print summary
print("Files with active GOTOs:\n")
for count in sorted(files_by_count.keys()):
    print(f"=== {count} GOTO(s) ({len(files_by_count[count])} files) ===")
    for f in sorted(files_by_count[count])[:10]:  # Show first 10
        print(f"  {f}")
    if len(files_by_count[count]) > 10:
        print(f"  ... and {len(files_by_count[count]) - 10} more")
    print()

total = sum(len(v) for v in files_by_count.values())
print(f"Total: {total} files with active GOTOs")
