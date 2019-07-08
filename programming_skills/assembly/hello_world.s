# Exit and send kernel a exit code
# Input: none
# Output: an exit status (use echo $? to check)
# Variables: %eax,%ebx
.section .data
.section .text
.globl _start
_start:
movl $1, %eax   # To exit this program (for kernel)
movl $0, %ebx   # Exit status
int $0x80
