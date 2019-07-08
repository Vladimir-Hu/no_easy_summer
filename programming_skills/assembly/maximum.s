 # Find the max value in the given numbers
 # Variables
 # %edi - Index of the given value
 # %ebx - Max number of each loop (also exit status)
 # %eax - Temp value
 .section .data
data_items:
 .long 3,69,34,222,45,75,55,36,66,33,22,11,66,78,98,41,0
 .section .text
 .globl _start
_start:
 movl $0,%edi
 movl data_items(,%edi,4),%eax
 movl %eax,%ebx
start_loop:
 cmpl $0,%eax
 je loop_exit
 incl %edi
 movl data_items(,%edi,4),%eax
 cmpl %ebx,%eax
 jle start_loop
 movl %eax,%ebx
 jmp start_loop
loop_exit:
 movl $1,%eax
 int $0x80
