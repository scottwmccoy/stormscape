"""available_scans must not lose today's volumes for a window ending today.

nexradaws lists day i as ``start + i days`` and skips it if that instant is in
the future, so 21 Oct 20Z -> 22 Oct 08Z listed on the morning of the 22nd
returned only the 21st. available_scans now asks one UTC day at a time.
"""
from datetime import datetime, timedelta, timezone

from stormscape import nexrad


class _Scan:
    def __init__(self, t):
        self.scan_time, self.filename = t, f"KRGX{t:%Y%m%d_%H%M%S}_V06"


class _BuggyConn:
    """Mimics nexradaws: per-day listing, days whose start+i is 'future' skipped."""
    def __init__(self, now):
        self.now, self.calls = now, []

    def get_avail_scans_in_range(self, start, end, radar):
        self.calls.append((start, end))
        out = []
        for i in range((end - start).days + 2):
            day = start + timedelta(days=i)
            if day > self.now:
                continue
            t = datetime(day.year, day.month, day.day, tzinfo=timezone.utc)
            while t.date() == day.date():
                if start <= t <= end:
                    out.append(_Scan(t))
                t += timedelta(minutes=30)
        return out


def test_window_ending_today_keeps_today(monkeypatch):
    now = datetime(2026, 10, 22, 15, tzinfo=timezone.utc)
    conn = _BuggyConn(now)
    monkeypatch.setattr(nexrad, "_conn", lambda: conn)
    start = datetime(2026, 10, 21, 20, tzinfo=timezone.utc)
    end = datetime(2026, 10, 22, 8, tzinfo=timezone.utc)
    # the bug, reproduced directly
    assert max(s.scan_time for s in conn.get_avail_scans_in_range(start, end, "KRGX")).day == 21
    scans = nexrad.available_scans("KRGX", start, end)
    times = [s.scan_time for s in scans]
    assert min(times) == start and max(times) == end
    assert len(times) == len(set(s.filename for s in scans)) == 25
