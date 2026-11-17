"""Checkpoints: the tutorial's 'line to check', as something that executes.

Every checkpoint scattered through the protocol exists because a command
produces plausible output even when the result is wrong. Prose cannot enforce
them; a human reading a 400-line LEaP log cannot either. Here each one is an
object that knows the question it asks, the evidence it read, and what it
concluded.

  PASS  the evidence says what it should
  WARN  proceed, but a human must have seen this
  FAIL  stop; anything downstream is meaningless

Nothing is ever skipped silently. A checkpoint whose evidence could not be
obtained is a FAIL, not an absence.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"
_RANK = {PASS: 0, WARN: 1, FAIL: 2}


@dataclass
class Check:
    name: str
    question: str
    result: str = PASS
    message: str = ""
    evidence: str = ""        # the actual line(s) read
    guidance: str = ""        # what to do about it
    step: str = ""

    def __str__(self):
        head = f"  [{self.result}] {self.name}: {self.message}"
        parts = [head]
        if self.evidence:
            ev = "\n".join("         | " + ln
                           for ln in self.evidence.strip().splitlines()[:6])
            parts.append(ev)
        if self.result != PASS and self.guidance:
            parts.append("         -> " + self.guidance)
        return "\n".join(parts)

    def to_dict(self):
        return asdict(self)


@dataclass
class CheckLog:
    checks: list = field(default_factory=list)

    def add(self, check: Check):
        self.checks.append(check)
        return check

    def record(self, name, question, ok, message, evidence="", guidance="",
               warn_only=False, step=""):
        res = PASS if ok else (WARN if warn_only else FAIL)
        return self.add(Check(name=name, question=question, result=res,
                              message=message, evidence=evidence,
                              guidance=guidance, step=step))

    @property
    def worst(self):
        return max((_RANK[c.result] for c in self.checks), default=0)

    @property
    def failed(self):
        return [c for c in self.checks if c.result == FAIL]

    @property
    def warned(self):
        return [c for c in self.checks if c.result == WARN]

    def ok(self):
        return not self.failed

    def summary(self):
        n = len(self.checks)
        return (f"{n} checkpoints: {n - len(self.failed) - len(self.warned)} "
                f"pass, {len(self.warned)} warn, {len(self.failed)} fail")

    def to_list(self):
        return [c.to_dict() for c in self.checks]
