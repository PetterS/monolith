#pragma once
#include <csignal>
#include <cstdio>
#include <functional>
#include <iostream>

#include <minimum/core/check.h>

#ifdef WIN32

namespace {

#define _WIN32_WINNT 0x0600
#define WIN32_LEAN_AND_MEAN
#include <assert.h>
#include <errno.h>
#include <process.h>
#include <stdio.h>
#include <windows.h>
#include <winnt.h>
#include <winternl.h>

struct SECTION_IMAGE_INFORMATION {
	PVOID EntryPoint;
	ULONG StackZeroBits;
	ULONG StackReserved;
	ULONG StackCommit;
	ULONG ImageSubsystem;
	WORD SubSystemVersionLow;
	WORD SubSystemVersionHigh;
	ULONG Unknown1;
	ULONG ImageCharacteristics;
	ULONG ImageMachineType;
	ULONG Unknown2[3];
};

struct RTL_USER_PROCESS_INFORMATION {
	ULONG Size;
	HANDLE Process;
	HANDLE Thread;
	CLIENT_ID ClientId;
	SECTION_IMAGE_INFORMATION ImageInformation;
};

#define RTL_CLONE_PROCESS_FLAGS_CREATE_SUSPENDED 0x00000001
#define RTL_CLONE_PROCESS_FLAGS_INHERIT_HANDLES 0x00000002
#define RTL_CLONE_PROCESS_FLAGS_NO_SYNCHRONIZE 0x00000004

#define RTL_CLONE_PARENT 0
#define RTL_CLONE_CHILD 297

typedef DWORD pid_t;

typedef NTSTATUS (*RtlCloneUserProcess_f)(
    ULONG ProcessFlags,
    PSECURITY_DESCRIPTOR ProcessSecurityDescriptor /* optional */,
    PSECURITY_DESCRIPTOR ThreadSecurityDescriptor /* optional */,
    HANDLE DebugPort /* optional */,
    RTL_USER_PROCESS_INFORMATION* ProcessInformation);

HANDLE child_process;

DWORD fork(void) {
	HMODULE mod = GetModuleHandle("ntdll.dll");
	if (!mod) {
		return -1;
	}

	RtlCloneUserProcess_f RtlCloneUserProcess =
	    (RtlCloneUserProcess_f)GetProcAddress(mod, "RtlCloneUserProcess");
	if (RtlCloneUserProcess == NULL) {
		return -1;
	}

	// fork.
	RTL_USER_PROCESS_INFORMATION process_info;
	memset(&process_info, 0, sizeof(process_info));
	process_info.Size = sizeof(process_info);
	NTSTATUS result = RtlCloneUserProcess(
	    RTL_CLONE_PROCESS_FLAGS_CREATE_SUSPENDED | RTL_CLONE_PROCESS_FLAGS_INHERIT_HANDLES,
	    NULL,
	    NULL,
	    NULL,
	    &process_info);

	if (result == RTL_CLONE_PARENT) {
		HANDLE me, ht, hcp = 0;
		DWORD pi, ti;
		me = GetCurrentProcess();
		pi = (DWORD)process_info.ClientId.UniqueProcess;
		ti = (DWORD)process_info.ClientId.UniqueThread;

		child_process = OpenProcess(PROCESS_ALL_ACCESS, FALSE, pi);
		minimum_core_assert(child_process);
		ht = OpenThread(THREAD_ALL_ACCESS, FALSE, ti);
		minimum_core_assert(ht);

		ResumeThread(ht);
		CloseHandle(ht);
		return pi;
	} else if (result == RTL_CLONE_CHILD) {
		return 0;
	}
	return -1;
}

int wait(int* status) {
	minimum_core_assert(WaitForSingleObject(child_process, INFINITE) == WAIT_OBJECT_0);
	DWORD exit_code = STILL_ACTIVE;
	minimum_core_assert(GetExitCodeProcess(child_process, &exit_code) != 0);
	minimum_core_assert(exit_code != STILL_ACTIVE);
	minimum_core_assert(CloseHandle(child_process) != 0);
	*status = exit_code;
	// std::cerr << "Child process exit code: 0x" << std::hex << exit_code << std::endl;
	return 0;
}

constexpr bool WIFEXITED(int status) { return status < 0xC0000000; }

constexpr int WEXITSTATUS(int status) { return status; }
}  // namespace

#else

#include <sys/wait.h>
#include <unistd.h>

#endif  // WIN32

// Performs a Catch2 check that "statement" kills/exits the process.
#define CHECK_DEATH(statement)                                                                 \
	do {                                                                                       \
		if (::internal::run_code([&]() { statement; }) == ::internal::SPECIAL_NODEATH_VALUE) { \
			FAIL("Expected death.");                                                           \
		} else {                                                                               \
			SUCCEED("Statement died.");                                                        \
		}                                                                                      \
	} while (false)

// Performs a Catch2 check that "statement" does not kill/exit the process.
#define CHECK_ALIVE(statement)                                                                 \
	do {                                                                                       \
		if (::internal::run_code([&]() { statement; }) == ::internal::SPECIAL_NODEATH_VALUE) { \
			SUCCEED("Statement did not kill process.");                                        \
		} else {                                                                               \
			FAIL("Statement killed process.");                                                 \
		}                                                                                      \
	} while (false)

#ifdef WIN32
// Windows does not support specific exit codes yet.
#define CHECK_EXIT_CODE(statement, code) \
	do {                                 \
	} while (false)
#else
// Performs a Catch2 check that "statement" exits the process with code.
#define CHECK_EXIT_CODE(statement, code)                                              \
	do {                                                                              \
		int actual_code = ::internal::run_code([&]() { statement; });                 \
		const bool process_died = actual_code == -1;                                  \
		CHECK_FALSE(process_died);                                                    \
		const bool process_exited = actual_code != ::internal::SPECIAL_NODEATH_VALUE; \
		CHECK(process_exited);                                                        \
		CHECK(actual_code == code);                                                   \
	} while (false)
#endif

namespace internal {
namespace {
constexpr int SPECIAL_NODEATH_VALUE = 112;

// Forks and calls f. Returns:
// * -1 if the process died, segfaulted etc.
// * The exit code if the process exited normally.
// * SPECIAL_NODEATH_VALUE if f ran until competion or threw
//   an exception.
int run_code(const std::function<void()>& f) {
	int result = fork();
	if (result == 0) {
		std::signal(SIGABRT, SIG_DFL);
		std::signal(SIGFPE, SIG_DFL);
		std::signal(SIGILL, SIG_DFL);
		std::signal(SIGINT, SIG_DFL);
		std::signal(SIGSEGV, SIG_DFL);
		std::signal(SIGTERM, SIG_DFL);

		std::fclose(stdout);
		std::fclose(stdin);

		try {
			f();
		} catch (...) {
			// Do nothing.
		}
#ifdef WIN32
		TerminateProcess(GetCurrentProcess(), SPECIAL_NODEATH_VALUE);
#else
		std::exit(SPECIAL_NODEATH_VALUE);
#endif
	} else if (result < 0) {
		FAIL("Fork failed.");
	} else {
		int status = 0;
		wait(&status);
		if (WIFEXITED(status)) {
			return WEXITSTATUS(status);
		}
	}
	return -1;
}
}  // namespace
}  // namespace internal
