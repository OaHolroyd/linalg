# compiler/linker
CC=/opt/homebrew/opt/llvm/bin/clang-20
LD=$(CC)

# libraries
INCLUDES=-I.
LDLIBS=-lm

# flags
WARNINGS=-Wall -Wextra -pedantic -Wno-unused-parameter -Wshadow \
         -Waggregate-return -Wbad-function-cast -Wcast-align -Wcast-qual \
         -Wfloat-equal -Wformat=2 -Wmissing-include-dirs \
         -Wnested-externs -Wpointer-arith -Wconversion -Wno-sign-conversion \
         -Wredundant-decls -Wsequence-point -Wstrict-prototypes -Wswitch -Wundef \
         -Wunused-but-set-parameter -Wwrite-strings -Wvla -Wno-gnu-zero-variadic-macro-arguments

# switch between debug and release mode
DEBUG ?=0
ifeq ($(DEBUG),0)
	CFLAGS=-O3
	LDFLAGS=
	PRE_CHECK=
	AT_CHECK=
	POST_CHECK=
else
	CFLAGS=-g3 -Og -fsanitize=address,undefined,leak -fno-omit-frame-pointer
	LDFLAGS=-g3 -fsanitize=address,undefined,leak

	# on macOS LSAN will detect leaks in the ObjC standard library so we have to suppress them
	PRE_CHECK=echo "leak:_fetchInitializingClassList\nleak:__Balloc\nleak:__parsefloat_buf" >> suppr.txt
	AT_CHECK=ASAN_OPTIONS=detect_leaks=1 LSAN_OPTIONS=suppressions=suppr.txt MallocNanoZone=0
	POST_CHECK=rm suppr.txt
endif

# directories
SRC_DIR=./src
OBJ_DIR=./obj
TEST_DIR=./tests

# files relating to the linear algebra code
SRC=$(wildcard $(SRC_DIR)/*.c)
OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.c=.o)))
DEPS=$(patsubst %.o,%.d,$(OBJ))

# files relating to test code
SRC_TEST=$(wildcard $(SRC_DIR)/tests/*.c)
OBJ_TEST=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_TEST:.c=.o)))
DEPS_TEST=$(patsubst %.o,%.d,$(OBJ_TEST))
TESTS=$(addprefix $(TEST_DIR)/, $(notdir $(SRC_TEST:.c=)))
RUN_TESTS=$(addprefix run_, $(notdir $(TESTS)))


# build the text editor executable
.PHONY: all
all: $(TESTS)

# return the project to its pre-build state
.PHONY: clean
clean:
	@printf "`tput bold``tput setaf 1`Cleaning`tput sgr0`\n"
	rm -rf $(OBJ_DIR) $(TEST_DIR)

# rebuild the project
.PHONY: rebuild
rebuild: | clean all

# run all tests
.PHONY: check
check: $(RUN_TESTS)

# link the linalg objects and the correct test objects to make the tests
$(TESTS): $(TEST_DIR)/% : $(OBJ_DIR)/%.o $(OBJ) | $(TEST_DIR)
	@printf "`tput bold``tput setaf 2`Linking %s`tput sgr0`\n" $@
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $< $(LDLIBS)

# compile object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@printf "`tput bold``tput setaf 6`Compiling %s`tput sgr0`\n" $@
	$(CC) $(CFLAGS) $(WARNINGS) $(INCLUDES) -MMD -MP -c -o $@ $<

# compile test objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/tests/%.c | $(OBJ_DIR)
	@printf "`tput bold``tput setaf 6`Compiling %s`tput sgr0`\n" $@
	$(CC) $(CFLAGS) $(WARNINGS) $(INCLUDES) -MMD -MP -c -o $@ $<

# include dependency information
-include $(DEPS)
-include $(DEPS_TEST)

# create directory for object files
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# create directory for test executables
$(TEST_DIR):
	mkdir -p $(TEST_DIR)

# run all test executables
.PHONY: $(RUN_TESTS)
$(RUN_TESTS): run_% : $(TESTS)
	@$(PRE_CHECK)
	@$(AT_CHECK) $(TEST_DIR)/$* \
	&& printf "`tput bold``tput setaf 2`PASSED %s`tput sgr0`\n" $* \
	|| printf "`tput bold``tput setaf 1`FAILED %s`tput sgr0`\n" $*
	@$(POST_CHECK)
