BUILD_DIR=/tmp/eggd_athena

.PHONY : all clean

all: create_tmp build_app clean_up

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
	@cp -r packages $(BUILD_DIR)/resources/home/dnanexus/

	# ensure output dir for app clear
	rm -f $(BUILD_DIR)/resources/home/dnanexus/athena/output/*

build_app:
	@echo "Building DNAnexus app"
	@dx build --app $(BUILD_DIR)

clean_up:
	@rm -r $(BUILD_DIR)