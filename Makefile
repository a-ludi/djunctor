build:
	dub build

test:
	dub test

format:
	find source/ -name '*.d' -type f -exec dfmt -i {} \;

clean:
	dub clean