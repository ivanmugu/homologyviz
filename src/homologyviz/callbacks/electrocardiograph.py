import time


class HeartBeatsParameters:
    """Parameters to monitor heart beats of the Dash app.

    The monitoring of the heart beats allows to stop the server when the app tab is closed
    in the browser.

    Attributes
    ----------
    last_heartbeat : dict
        A dictionary storing the timestamp of the last heartbeat and a counter.
    timeout_seconds : int
        The number of seconds before a timeout occurs if no heartbeat is received.
    heartbeat_monitor_started : bool
        Whether the heartbeat monitor has been started
    """

    def __init__(
        self,
        last_heartbeat: dict | None = None,
        timeout_seconds: int = 5,
        heartbeat_monitor_started: bool = False,
    ) -> None:
        """Initialize HeartBeatsParameters

        Parameters
        ----------
        last_heartbeat : dict, optional
            Initial dictionary storing the timestamp and counter. Defaults to current time.
        timeout_seconds : int, optional
            Timeout duration in seconds. Default is 5 seconds.
        heartbeat_monitor_started : bool, optional
            Whether the monitor is started. Default is False.
        """
        self.last_heartbeat = (
            last_heartbeat
            if last_heartbeat is not None
            else {"timestamp": time.time(), "counter": 0}
        )
        self.timeout_seconds = timeout_seconds
        self.heartbeat_monitor_started = heartbeat_monitor_started
