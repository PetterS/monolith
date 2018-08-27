
#ifdef EMSCRIPTEN
// This fucntion from libc ends up not being included by Emscripted.
// Having it here removes the need for manual editing of the resulting
// JS files.
extern "C" {
int __cxa_thread_atexit(void (*)(char*), char*, char*) { return 0; }
}
#endif
