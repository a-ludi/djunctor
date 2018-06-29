build:
	dub build

test:
	dub test

format:
	find source/ -name '*.d' -type f -exec dfmt -i {} \; -exec sed -i 's/\s\+$$//' {} \;

clean:
	dub clean

todos:
	@grep --color --exclude='*.tar'{,.gz,.xz} -noPR '\b(TODO|FIXME)[\s]*?[:\s][\s]*(?P<todo>.*)$$' ./source ./tests
