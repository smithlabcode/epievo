all_subdirs=prog ${SMITHLAB_CPP} ${TREETOOL}
app_subdirs=prog

all:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
	make -C $${i} SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) test; \
	done;

clean:
	@for i in $(all_subdirs); do \
	make -C $${i} clean; \
	done;

.PHONY: test
