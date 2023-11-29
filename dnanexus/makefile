BUILD_DIR=/tmp/eggd_athena

.PHONY : all clean

all: create_tmp build_packages build_app clean_up

create_tmp:
	@echo "Copying files"
	@[ -d $(BUILD_DIR) ] && rm -r $(BUILD_DIR)/* || mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/resources/home/dnanexus
	@mkdir -p $(BUILD_DIR)/src

	# copy in app code
	@cp -r dxapp.json readme.md $(BUILD_DIR)/
	@cp -r src/eggd_athena.sh $(BUILD_DIR)/src
	@cp -r resources/bedtools.static.binary.gz $(BUILD_DIR)/resources/home/dnanexus/
	@cp -r ../athena $(BUILD_DIR)/resources/home/dnanexus/

	# get required Python packages
	mkdir -p $(BUILD_DIR)/resources/home/dnanexus/packages
	pip download -r ../requirements.txt -d $(BUILD_DIR)/resources/home/dnanexus/packages

build_packages:
	# ensure any Python source packages built as .whl files instead of source
	# code for faster installing when running the app
	for pkg in $(BUILD_DIR)/resources/home/dnanexus/packages/*.tar.gz; do\
		mkdir $(BUILD_DIR)/resources/home/dnanexus/packages/build;\
		tar xf $${pkg} -C $(BUILD_DIR)/resources/home/dnanexus/packages/build/;\
		pip wheel --no-deps -w $(BUILD_DIR)/resources/home/dnanexus/packages/ $${pkg};\
		rm $${pkg};\
		rm -r $(BUILD_DIR)/resources/home/dnanexus/packages/build;\
	done

build_app:
	@echo "Building DNAnexus app"
	@dx build --app $(BUILD_DIR)

clean_up:
	@rm -r $(BUILD_DIR)
