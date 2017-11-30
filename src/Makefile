all_subdirs=prog common ${SMITHLAB_CPP}
app_subdirs=prog

all:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) test; \
	done;

clean:
	@for i in $(all_subdirs); do \
	make -C $${i} clean; \
	done;

.PHONY: test
