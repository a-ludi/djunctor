.PHONY: build test format stylecheck unusedimports clean todos

build:
	dub build

test:
	dub test

format:
	find source/ -name '*.d' -type f -exec dfmt -i {} \; -exec sed -i 's/\s\+$$//' {} \;

stylecheck:
	@dub run dscanner -- --styleCheck | less

unusedimports:
	@./.unused-imports.sh $$(find source/ -name '*.d' -type f)

clean:
	dub clean

todos:
	@grep --color --exclude='*.tar'{,.gz,.xz} -noPR '\b(TODO|FIXME)[\s]*?[:\s][\s]*(?P<todo>.*)$$' ./source ./tests
