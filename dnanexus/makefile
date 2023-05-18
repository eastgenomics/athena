BUILD_DIR=/tmp/eggd_athena

.PHONY : all clean

all: create_tmp build_app clean_up

create_tmp:
	@echo "Copying files" 
	@[ -d $(BUILD_DIR) ] && rm -r $(BUILD_DIR)/* || mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/resources/home/dnanexus

	@cp -r dxapp.json eggd_athena.sh readme.md $(BUILD_DIR)
	@cp -r ../athena $(BUILD_DIR)/resources/home/dnanexus

build_app:
	@echo "Building DNAnexus app"
	@dx build --app $(BUILD_DIR)

clean_up:
	@rm -r $(BUILD_DIR)
