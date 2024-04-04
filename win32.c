#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#define OS_READ     GENERIC_READ
#define OS_WRITE    GENERIC_WRITE
#define OS_RW       (GENERIC_READ | GENERIC_WRITE)
#define OS_SEEK_BEG FILE_BEGIN
#define OS_SEEK_CUR FILE_CURRENT
#define OS_SEEK_END FILE_END

typedef HANDLE os_file;

static os_file
os_open(s8 path, int flags)
{
	os_file f;
	f = CreateFileA((char *)path.data, flags, 0, NULL,
	                OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (f == INVALID_HANDLE_VALUE)
		die("rip can't open output file: %s\n", path.data);
	return f;
}

static void
os_write(os_file f, s8 s)
{
	if (WriteFile(f, s.data, s.len, 0, 0) == 0)
		die("can't write to file\n");
}

static void
os_close(os_file f)
{
	CloseHandle(f);
}

static void
os_seek(os_file f, size off, int whence)
{
	SetFilePointer(f, off, 0, whence);
}

static u32
os_get_core_count(void)
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
}

static void
os_pause(void)
{
	fputs("Press any key to close...", stdout);
	fgetc(stdin);
}
