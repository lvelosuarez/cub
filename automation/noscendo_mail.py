#!/usr/bin/env python3
"""noscendo_mail.py <subject> <body> [recipient] — send via CHU relay."""
import smtplib, sys
from email.message import EmailMessage
SMTP_HOST, SMTP_PORT = "smtp.chu-brest.fr", 25
SENDER = "cub-pipeline@chu-brest.fr"
def main():
    if len(sys.argv) < 3:
        print("usage: noscendo_mail.py <subject> <body> [recipient]", file=sys.stderr); return 2
    msg = EmailMessage()
    msg["From"] = SENDER
    msg["To"] = sys.argv[3] if len(sys.argv) > 3 else "lourdesvelo@gmail.com"
    msg["Subject"] = sys.argv[1]
    msg.set_content(sys.argv[2])
    with smtplib.SMTP(SMTP_HOST, SMTP_PORT, timeout=20) as s:
        s.send_message(msg)
    print("sent ok"); return 0
if __name__ == "__main__":
    sys.exit(main())
