BITS 64
; /**
;  * @file bignum_jacobi.asm
;  * @brief System V AMD64 implementation of the public Jacobi-symbol API.
;  * @details
;  * This file implements bignum_jacobi() and its private binary-reduction
;  * helper. The algorithm accepts non-negative little-endian bignum_t values,
;  * removes powers of two from the numerator, applies the supplementary law
;  * for 2 and quadratic reciprocity, and reduces by binary long division.
;  *
;  * C/ASM boundary (Linux x86-64, System V AMD64 ABI):
;  *   RDI = const bignum_t *a      (borrowed input numerator)
;  *   RSI = const bignum_t *n      (borrowed input positive odd modulus)
;  *   RDX = int *symbol            (caller-allocated output)
;  *   EAX = bignum_jacobi_status_t return value
;  * bignum_t is represented as 32 little-endian uint64_t words at byte offset
;  * 0, followed by uint64_t len at byte offset 256. len is normalized and is
;  * in the inclusive range 0..BIGNUM_CAPACITY. The output is written only on
;  * successful validation and receives -1, 0, or 1.
;  *
;  * Register contract:
;  *   Preserved: RBX, RBP, R12, R13. These are saved/restored by the public
;  *   entry point and by jacobi_reduce when it is called internally.
;  *   Caller-saved/clobbered: RAX, RCX, RDX, RSI, RDI, R8, R9, R10, R11.
;  *   RDX is used as the output pointer after entry; no condition flags are
;  *   part of the public interface. XMM/YMM state is not used or modified.
;  *
;  * Stack frame:
;  *   The entry RSP is 8 mod 16. Four pushes save RBX/RBP/R12/R13 and the
;  *   public function reserves 1096 bytes; consequently RSP is 0 mod 16
;  *   immediately before every CALL, as required by the ABI.
;  *   After allocation, the public frame uses:
;  *     [rsp+0x000 .. rsp+0x107]  local copy of a (33 qwords)
;  *     [rsp+0x110 .. rsp+0x217]  mutable numerator state
;  *     [rsp+0x220 .. rsp+0x327]  mutable modulus state
;  *     [rsp+0x330 .. rsp+0x437]  swap temporary state
;  *   The final 16 bytes of the reserved area maintain frame separation; saved
;  *   callee-saved registers are above the allocated frame.
;  *
;  * Error semantics: NULL pointers return BIGNUM_JACOBI_ERROR_NULL_ARG;
;  * zero/even modulus returns BIGNUM_JACOBI_ERROR_MODULUS; all error paths
;  * leave the caller-owned output unchanged. The implementation allocates no
;  * heap memory and is reentrant/thread-safe for independently owned inputs.
;  *
;  * The private jacobi_reduce helper receives RDI=x and RSI=m, where both are
;  * writable local bignum_t states. It preserves the same callee-saved
;  * registers, uses a stack-local 33-word remainder, and writes the reduced
;  * value back to x. It is not a public ABI symbol.
;  */

BITS 64
DEFAULT REL

jacobi_reduce:
push	r13
mov	r9, rsi
mov	r8, rdi
mov	ecx, 33
push	r12
xor	eax, eax
push	rbp
push	rbx
sub	rsp, 160
lea	rsi, [rsp-120]
mov	rdi, rsi
rep stosq
mov	rdi, qword [r8+256]
test	rdi, rdi
je	L1
sal	rdi, 6
je	L3
mov	r11, qword [r9+256]
lea	r10, [rdi-1]
xor	edi, edi
L21:
mov	rax, r10
mov	ecx, r10d
shr	rax, 6
mov	rax, qword [r8+rax*8]
shr	rax, cl
mov	rcx, rax
and	ecx, 1
test	rdi, rdi
je	L4
mov	rdx, rsi
lea	rbp, [rsi+rdi*8]
L5:
mov	rax, qword [rdx]
mov	rbx, rcx
add	rdx, 8
mov	rcx, rax
add	rax, rax
or	rax, rbx
shr	rcx, 63
mov	qword [rdx-8], rax
cmp	rdx, rbp
jne	L5
L4:
cmp	rdi, 32
je	L6
test	rcx, rcx
je	L6
mov	qword [rsp+rdi*8-120], 1
lea	rdx, [rdi+1]
cmp	rdx, r11
je	L7
mov	rdi, rdx
cmp	rdx, r11
jb	L8
L9:
test	r11, r11
je	L24
cmp	rdx, r11
mov	r12, r11
cmovbe	r12, rdx
xor	edi, edi
xor	eax, eax
L15:
mov	rbp, qword [rsi+rax*8]
mov	rcx, qword [r9+rax*8]
mov	rbx, rbp
sub	rbx, rcx
mov	r13, rbx
sub	r13, rdi
cmp	rbp, rcx
setb	cl
cmp	rbx, rdi
mov	qword [rsi+rax*8], r13
setb	dil
add	rax, 1
or	ecx, edi
movzx	edi, cl
cmp	rax, r12
jb	L15
cmp	rax, rdx
jnb	L19
L17:
mov	rcx, qword [rsi+rax*8]
cmp	rcx, rdi
setb	bl
sub	rcx, rdi
mov	qword [rsi+rax*8], rcx
add	rax, 1
movzx	edi, bl
cmp	rax, rdx
jb	L17
mov	rdi, rdx
sub	rdx, 1
cmp	qword [rsi+rdx*8], 0
jne	L8
L20:
test	rdx, rdx
je	L13
L19:
mov	rdi, rdx
sub	rdx, 1
cmp	qword [rsi+rdx*8], 0
je	L20
L8:
sub	r10, 1
jnb	L21
L3:
mov	qword [rsp+136], rdi
mov	ecx, 33
mov	rdi, r8
rep movsq
L1:
add	rsp, 160
pop	rbx
pop	rbp
pop	r12
pop	r13
ret
L6:
cmp	rdi, r11
je	L7
mov	rdx, rdi
jb	L8
L10:
test	rdx, rdx
jne	L9
L13:
xor	edi, edi
sub	r10, 1
jnb	L21
jmp	L3
L7:
mov	rax, r11
jmp	L11
L12:
sub	rax, 1
mov	rbx, qword [r9+rax*8]
cmp	qword [rsi+rax*8], rbx
jne	L54
L11:
test	rax, rax
jne	L12
L51:
mov	rdx, r11
jmp	L10
L24:
xor	edi, edi
xor	eax, eax
jmp	L17
L54:
mov	rdi, r11
jnb	L51
sub	r10, 1
jnb	L21
jmp	L3
global bignum_jacobi
bignum_jacobi:
push	r13
push	r12
push	rbp
push	rbx
mov	rbx, rdx
sub	rsp, 1096
test	rsi, rsi
sete	al
test	rdx, rdx
sete	dl
or	eax, edx
test	rdi, rdi
sete	dl
or	al, dl
jne	L73
mov	r8, rsi
mov	rdx, rsp
mov	r9, rdi
mov	ecx, 33
mov	rax, qword [r8+256]
mov	rdi, rdx
rep movsq
test	rax, rax
jne	L57
jmp	L78
L58:
test	rax, rax
je	L78
L57:
mov	r8, rax
sub	rax, 1
cmp	qword [rdx+rax*8], 0
je	L58
test	byte [rsp], 1
je	L78
lea	r13, [rsp+544]
lea	rbp, [rsp+272]
mov	rsi, r9
mov	qword [rsp+256], r8
lea	rdi, [rsp+272]
mov	ecx, 33
mov	r12d, 1
rep movsq
mov	rsi, rdx
lea	rdi, [rsp+544]
mov	ecx, 33
rep movsq
mov	rsi, r13
mov	rdi, rbp
call	jacobi_reduce
mov	rax, qword [rsp+528]
test	rax, rax
je	L61
L60:
mov	r10, qword [rsp+544]
xor	r9d, r9d
mov	rsi, r10
and	esi, 7
sub	rsi, 3
and	rsi, -3
test	rax, rax
je	L61
L101:
mov	rdx, qword [rsp+272]
test	dl, 1
jne	L99
lea	rdx, [rbp+rax*8+0]
xor	r8d, r8d
L62:
mov	rdi, qword [rdx-8]
sub	rdx, 8
mov	rcx, rdi
sal	rdi, 63
shr rcx, 1
or	rcx, r8
mov	r8, rdi
mov	qword [rdx], rcx
cmp	rbp, rdx
jne	L62
mov	rdx, rax
sub	rax, 1
cmp	qword [rbp+rax*8+0], 0
jne	L100
L65:
mov	r9d, 1
test	rax, rax
je	L64
mov	rdx, rax
sub	rax, 1
cmp	qword [rbp+rax*8+0], 0
je	L65
L100:
mov	rax, rdx
L64:
mov	edx, r12d
neg	r12d
test	rsi, rsi
cmovne	r12d, edx
test	rax, rax
jne	L101
L61:
cmp	qword [rsp+800], 1
je	L102
L77:
xor	r12d, r12d
L71:
mov	dword [rbx], r12d
xor	eax, eax
L55:
add	rsp, 1096
pop	rbx
pop	rbp
pop	r12
pop	r13
ret
L99:
test	r9b, r9b
je	L68
mov	qword [rsp+528], rax
L68:
not	rdx
and	edx, 3
jne	L69
mov	eax, r12d
not	r10
neg	eax
and	r10d, 3
cmove	r12d, eax
L69:
mov	rsi, rbp
lea	rdi, [rsp+816]
mov	ecx, 33
rep movsq
lea	rbp, [rsp+272]
mov	rsi, r13
mov	ecx, 33
lea	r13, [rsp+544]
lea	rdi, [rsp+272]
rep movsq
lea	rdi, [rsp+544]
lea	rsi, [rsp+816]
mov	ecx, 33
rep movsq
mov	rsi, r13
mov	rdi, rbp
call	jacobi_reduce
mov	rax, qword [rsp+528]
test	rax, rax
jne	L60
cmp	qword [rsp+800], 1
jne	L77
L102:
xor	eax, eax
cmp	qword [rsp+544], 1
cmovne	r12d, eax
jmp	L71
L78:
add	rsp, 1096
mov	eax, -2
pop	rbx
pop	rbp
pop	r12
pop	r13
ret
L73:
mov	eax, -1
jmp	L55
