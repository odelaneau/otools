projects = liftover mendel swapalleles

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm -f static_bins/*
	rm -f docker/resources/*
	rm -f docker/shapeit5*.tar.gz

