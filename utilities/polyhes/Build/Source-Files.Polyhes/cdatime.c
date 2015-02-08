/* this is a c program to call ctime and return the ascii date string */
#include <sys/types.h>

void cdate(tmpstr,sec_time)
char *tmpstr;
long int *sec_time;

{

   char *ctime();
   strcpy(tmpstr,ctime(sec_time));
}
/* end of c routine */

