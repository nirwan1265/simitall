"""Minimal OBO parser shared by the lexicon build scripts."""
import gzip, re

SYN = re.compile(r'^synonym: "((?:[^"\\]|\\.)*)" (\w+)(?: ([\w:]+))?')


def parse(path):
    op = gzip.open if path.endswith(".gz") else open
    terms, cur = {}, None
    with op(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("["):
                cur = {"id": None, "name": None, "syn": [], "is_a": [], "obsolete": False} if line == "[Term]" else None
                continue
            if cur is None or not line:
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:]
                terms[cur["id"]] = cur
            elif line.startswith("name: "):
                cur["name"] = line[6:]
            elif line.startswith("synonym: "):
                m = SYN.match(line)
                if m:
                    cur["syn"].append((m.group(1).replace('\\"', '"'), m.group(2), m.group(3) or ""))
            elif line.startswith("is_a: "):
                cur["is_a"].append(line[6:].split(" ")[0])
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    return terms
