global distance8
global maxdist

distance8:
        push ebx
        push edi
        push esi
        push edx
        push ecx
        push eax
        
        mov ecx,[dword ESP+24+28]
        and cl,127
        xor edi,edi
fillmax:
        mov [maxdist+edi], cl
        inc edi
        cmp edi, 16
        jne fillmax

        xor eax,eax ; position in the array
        xorpd xmm5,xmm5 ; count of changes
again:
        movdqa xmm4, [b255] ; best cluster
        movdqa xmm7, [b127] ; best distance
        xorpd xmm6, xmm6 ; current cluster number (16-fold)
        xor edx,edx
myloop2: 
        xorpd xmm1, xmm1 ; max values
        mov edi,[dword ESP+24+8] ; clusters
        mov ecx,[dword ESP+24+12] ; num of dimensions
        mov edi, [edi+edx*4] ; cluster data (16-fold)
        mov esi,[dword ESP+24+4] ; points
myloop:
        ; load first dimension of 16 points
        mov ebx, [esi]
        add ebx, eax
        add esi, 4
        movdqa xmm0, [ebx]
        ; compute absolute difference of current dimension
        psubb xmm0, [edi]
        add edi, 16

        ; replacement for
        ; psign xmm0, xmm0 
        ; or
        ; pabsb xmm0,xmm0
        ; which nasm doesn't know
        xorpd xmm2,xmm2
        pcmpgtb xmm2, xmm0
        xorpd xmm0, xmm2
        psubb xmm0, xmm2

        pmaxub xmm1, xmm0
        dec ecx
        jnz myloop

        ; now, xmm1 contains the maximum absolute distances

        ; remove any differences which are larger than the current maximum
        movdqa xmm3, xmm1
        pcmpgtb xmm3, [maxdist]
        pxor xmm1, [b128]
        por xmm1, xmm3 ; set all distances to 127 which are larger than maxdist
        pxor xmm1, [b128]

        ; figure out whether the distance is better than the currently best cluster
        movdqa xmm2, xmm7
        pcmpgtb xmm2, xmm1 ; now xmm2 contains the byte mask of all points where this center is better (SMALLER!)
       
        ; move current cluster (xmm3) into xmm4 where xmm2 says it's better
        movdqa xmm3, xmm6
        pand xmm3, xmm2
        pxor xmm2, [b255] ; there doesn't seem to be a 128-bit neg in sse2
        pand xmm4, xmm2
        por xmm4, xmm3

        pminub xmm7, xmm1
        paddb xmm6, [one]
        
        inc edx
        cmp edx, [dword ESP+24+16] ; compare with num of clusters
        jne myloop2

        mov edi,[dword ESP+24+20] ; dest

        movdqa xmm0,xmm4
        pcmpeqb xmm0,[edi] ; find all clusters which didn't change
        pxor xmm0, [one]
        pand xmm0, [one] ; 01 for each cluster which changed
        paddusb xmm5, xmm0 ; calculate count of changes

        movaps [edi], xmm4
        add edi, 16
        add eax, 16
        mov [dword ESP+24+20],edi
        mov ecx, [dword ESP+24+24]
        dec ecx
        mov [dword ESP+24+24],ecx
        jnz again

        ;movdqu [tmp],xmm5
        movdqa [tmp],xmm5 ; tmp is aligned

        pop eax
        pop ecx
        pop edx
        pop esi

        xor ebx,ebx
        xor edi,edi
        xor eax,eax
addup:
        mov bl,[tmp+edi]
        add eax,ebx
        inc edi
        cmp edi, 16
        jne addup

        pop edi
        pop ebx
        ret

align 16
b127: db 127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127
b128: db 128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128
one:  db 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1
b255: db 255,255,255,255, 255,255,255,255, 255,255,255,255, 255,255,255,255

section .data align=16
align 16

tmp:  db 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
maxdist: db 8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8

;psadbw: sum of total absolute differences

