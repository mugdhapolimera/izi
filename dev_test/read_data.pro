PRO read_data, H, filename,COLUMNS=cols,ROWS=rows 
OPENR,1,filename 
IF N_ELEMENTS(cols) LE 0 THEN cols=1 ;Default value for cols 
IF N_ELEMENTS(rows) LE 0 THEN rows=1000    ;Default value for rows 
H=FLTARR(cols,rows) ;A big array to hold the data 
S=FLTARR(cols)      ;A small array to read a line 
ON_IOERROR,ers     ;Jump to statement ers when I/O error is detected 
n=0 ; Create a counter 
WHILE n LT rows DO BEGIN 
    READF,1,S    ;Read a line of data 
    H[*,n]=S     ;Store it in H 
    n=n+1        ;Increment the counter 
ENDWHILE          ;End of while loop 
ers: CLOSE,1         ;Jump to this statement when an end of file is detected 
H=H[*,0:n-1] 
END 