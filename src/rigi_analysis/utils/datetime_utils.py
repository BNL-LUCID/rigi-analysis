import datetime


def datetime_now() -> str:
    """Return current datetime as a formatted string."""
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
