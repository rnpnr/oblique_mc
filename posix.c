#include <fcntl.h>
#include <unistd.h>

#define OS_READ     O_RDONLY
#define OS_WRITE    (O_WRONLY | O_CREAT | O_TRUNC)
#define OS_RW       O_RDWR
#define OS_SEEK_CUR SEEK_CUR
#define OS_SEEK_END SEEK_END
#define OS_SEEK_SET SEEK_SET

typedef int os_file;

static os_file
os_open(s8 path, int flags)
{
	os_file f;
	f = open((char *)path.data, flags, 0600);
	if (f == -1)
		die("rip can't open output file: %s\n", path.data);
	return f;
}

static void
os_write(os_file f, s8 s)
{
	if (write(f, s.data, s.len) == -1)
		die("can't write to file\n");
}

static void
os_close(os_file f)
{
	close(f);
}

static void
os_seek(os_file f, size off, int whence)
{
	lseek(f, off, whence);
}

static u32
os_get_core_count(void)
{
	return sysconf(_SC_NPROCESSORS_ONLN);
}
