/* eliminate the '_' suffix */
/* cgetmem_(buf,len,pt,bpt) */
cgetmem(buf,len,pt,bpt)
int len[],buf[];
int *pt[],*bpt[];
{
 int *malloc();
 int *ptr, *ib;
 int ln;
 ln=len[0];
 ptr = malloc(ln);
 pt[0] = ptr;
 ib = &buf[0];
 bpt[0] = ib;
 return;
}
