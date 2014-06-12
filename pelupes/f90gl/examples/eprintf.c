/* On my SunOS workstation, the loader gives an unresolved reference
   to __eprintf.  This is just to resolve that reference.  Hopefully,
   you won't need it. */

void __eprintf(void){}
