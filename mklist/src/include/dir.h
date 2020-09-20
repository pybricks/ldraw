#ifndef DIR_H
#define DIR_H

#  include <dirent.h>
#  include <sys/stat.h>
#  include <unistd.h>

struct ffblk 
{
  char *ff_name;
};

#if defined(MAC)
  char separator = ':';
# define SEPSTR ":"
#else
  char separator = '/';
# define SEPSTR "/"
#endif

DIR *dirp;
struct dirent *dir;
char directory[300];

int findnext(struct ffblk *ffb)
{
  struct stat		statbuf;
  char	filename[300];

  for (dir = readdir(dirp); dir; dir = readdir(dirp))
  {
    ffb->ff_name = dir->d_name;
    strcpy(filename, directory);
    strcat(filename, SEPSTR);
    strcat(filename, dir->d_name);
    stat(filename,&statbuf);
    if ((statbuf.st_mode & S_IFDIR) == 0) // Skip directories
      return 0; // not done
  }
  return 1; // done
}

int findfirst(char *path, struct ffblk *ffb, int first)
{
  int i;
  char *ptr;

  /* Localize Directory Separators */
  for(i=0; i<strlen(path); i++) 
  {
    if ((path[i] == '/') || (path[i] == '\\'))
    {
      path[i] = separator; 
    }
    else 
      path[i] = tolower(path[i]);
  }

  strcpy(directory,path);
  /* Remove any trailing separators (and anything after them) */
  if ( (ptr = strrchr(directory, separator)) )
    *ptr = 0;

  dirp = opendir(directory);
  if (!dirp)
    return 1; // done

  return findnext(ffb);
}

#endif /* DIR_H */

