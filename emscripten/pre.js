// Begin pre.js.

// Send a message to the main thread.
var send_custom_message = function(command, data) {
  data = data || {};
  data["target"] = "custom";
  data["command"] = command;
  postMessage(data);
};

// Used to communicate with the main thread.
var custom_message_handler = function(event) {
  data = event.data.userData;
  if (data.command == "callMain") {
    Module["callMain"](data.arguments);
    // When we are done executing, post a completion message.
    postMessage({ target: "custom", command: "callMainComplete" });
  } else if (data.command == "FS_unlink") {
    Module["FS_unlink"](data.filename);
    console.warn("Removed file " + data.filename);
  } else if (data.command == "FS_createDataFile") {
    Module["FS_createDataFile"](
      data.dir,
      data.filename,
      data.content,
      true,
      true
    );
    console.warn("Created file " + data.filename + " in " + data.dir);
  } else if (data.command == "custom-init") {
    if (ENVIRONMENT_IS_NODE) {
      FS.mkdir("/db");
      FS.mount(NODEFS, { root: "." }, "/db");
    } else {
      FS.mkdir("/db");
      FS.mount(IDBFS, {}, "/db");
    }
    FS.syncfs(true, function(err) {
      if (err) {
        console.error(err);
      } else {
        console.warn("Mounted IndexedDB as a file system at /db.");
        send_custom_message("custom-init-done");
      }
    });
  } else if (data.command == "FS_sync") {
    FS.syncfs(false, function(err) {
      if (err) {
        console.error(err);
      } else {
        console.warn("Synced /db back to IndexedDB.");
      }
    });
  } else if (data.command == "ls") {
    FS.nameTable.forEach(function(node) {
      if (node && typeof node.parent !== "string") {
        var name = FS.getPath(node);
        send_custom_message("file_name", { name: name });
      } else {
        console.warn("Not including ", node);
      }
    });
  } else if (data.command == "cat") {
    var content = FS.readFile(data.name, { encoding: "utf8" });
    send_custom_message("cat", { name: data.name, content: content });
  } else if (data.command == "write-file") {
    var content = FS.writeFile(data.name, data.content, { encoding: "utf8" });
  } else if (data.command == "delete-file") {
    FS.unlink(data.name);
  } else {
    console.error("Unknown custom command: ", data);
  }
};

// Create a module that does not run the program by default.
var Module = {
  noInitialRun: true,
  noExitRuntime: true,
  onCustomMessage: custom_message_handler
};

setTimeout(function() {
  console.warn("Web worker script loaded.");
}, 0);

// End pre.js.
