import json
import os
import signal
import threading
import time

from flask import Flask, Response, jsonify, request

from homologyviz.callbacks.electrocardiograph import HeartBeatsParameters


def register_server_lifecycle(
    server: Flask,
    heartbeat_parameters: HeartBeatsParameters,
) -> None:
    """Register heartbeat monitoring and shutdown routes."""

    @server.route("/heartbeat", methods=["POST"])
    def heartbeat() -> tuple[Response, int]:
        """
        Receive heartbeat pings from the frontend to monitor whether the app tab is open.

        This route is periodically called by the frontend to indicate that the app is
        still active. It parses the POST payload (JSON or raw) and updates the internal
        heartbeat counter and timestamp. If no data is received, it returns a warning.

        Returns
        -------
        tuple
            A Flask response with a JSON payload indicating success or failure,
            and an HTTP status code (200 or 500).
        """
        try:
            data = None

            # Attempt to parse the JSON payload
            if request.is_json:
                data = request.get_json()
            elif request.data:
                data = json.loads(request.data.decode("utf-8"))

            # Handle cases where no data is received
            if not data:
                print("Warning: No data received in the heartbeat request.", flush=True)
                return jsonify(success=False, message="No data received"), 200

            counter = data.get("counter", 0)
            heartbeat_parameters.last_heartbeat["timestamp"] = time.time()
            heartbeat_parameters.last_heartbeat["counter"] = counter

            return jsonify(success=True), 200
        except Exception as e:
            print(f"Error in /heartbeat route: {e}", flush=True)
            return jsonify(success=False, error=str(e)), 500

    @server.route("/shutdown", methods=["POST"])
    def shutdown_server() -> tuple[str, int]:
        """
        Shut down the Dash server when triggered.

        This endpoint is called by `monitor_heartbeats` when the app tab is closed
        and no heartbeats are received for a prolonged period. It sends a SIGINT
        signal to terminate the current process.

        Returns
        -------
        tuple
            A string message and HTTP status code 200 indicating shutdown.
        """
        os.kill(os.getpid(), signal.SIGINT)  # Send a signal to terminate the process
        print("Server shutting down...")
        return "Server shutting down...", 200

    def monitor_heartbeats() -> None:
        """
        Continuously monitor heartbeat timestamps to detect tab closure and shut down the
        server.

        This function runs in a background thread and checks whether the most recent
        heartbeat has timed out (based on `heartbeat_parameters.timeout_seconds`).
        If no new heartbeat is detected for a set period, and the heartbeat counter
        remains unchanged, the server is shut down gracefully.
        """
        counter = 0
        while True:
            now = time.time()
            elapsed_time = now - heartbeat_parameters.last_heartbeat["timestamp"]
            counter += 1
            # If timeout occurs, shut down the server
            if elapsed_time > heartbeat_parameters.timeout_seconds:
                print("Timeout: No heartbeats. Checking if counter has stopped...")
                # Check if the counter has stopped increasing
                initial_counter = heartbeat_parameters.last_heartbeat["counter"]
                time.sleep(5)  # Wait to see if the counter increases
                if heartbeat_parameters.last_heartbeat["counter"] == initial_counter:
                    shutdown_server()
            time.sleep(1)  # Regular monitoring interval

    if not heartbeat_parameters.heartbeat_monitor_started:
        heartbeat_parameters.heartbeat_monitor_started = True
        print("Initiating heartbeat_monitor_started")
        # Start the monitoring thread
        threading.Thread(
            target=monitor_heartbeats,
            daemon=True,
        ).start()
