import re

from detect_secrets.plugins.base import RegexBasedDetector

class OCCRegexDetector(RegexBasedDetector):
    """Scans for common secrets in OCC."""
    secret_type = 'OCC REGEX'


    denylist = (
        re.compile(
            r"([0-9a-zA-Z\/:=,&?_]*maps\.googleapis\.com[0-9a-zA-Z\/:=,&?_|${}]*key=[0-9a-zA-Z]+)",
            re.IGNORECASE,
        ),
    )
