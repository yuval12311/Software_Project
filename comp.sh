gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c vector.c
gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c spkmeans.c
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.o vector.o -lm -o spkmeans
