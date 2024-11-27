let heartbeatCounter = 0;

function sendHeartbeat() {
    heartbeatCounter += 1; // Increment the counter
    console.log("Heartbeat sent with counter:", heartbeatCounter);
    const data = JSON.stringify({ counter: heartbeatCounter });
    const blob = new Blob([data], { type: 'application/json' });
    navigator.sendBeacon('/heartbeat', blob);
}

// Send a heartbeat every 5 seconds
setInterval(sendHeartbeat, 5000);

// On tab unload, send a final counter value (optional)
window.addEventListener("beforeunload", function () {
    navigator.sendBeacon('/heartbeat', JSON.stringify({ counter: heartbeatCounter }));
});
