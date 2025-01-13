BUILD_DIR=/tmp/eggd_athena

.PHONY : all clean

all: create_tmp get_python_packages get_bedtools_binary build_app clean_up

create_tmp:
	@printf "\n\nCopying files for build"
	@[ -d $(BUILD_DIR) ] && rm -r $(BUILD_DIR)/* || mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/resources/home/dnanexus/packages
	@mkdir -p $(BUILD_DIR)/resources/home/dnanexus/athena
	@mkdir -p $(BUILD_DIR)/src

	# copy in app code
	@cp -r dxapp.json readme.md $(BUILD_DIR)/
	@cp -r src/eggd_athena.sh $(BUILD_DIR)/src
	@cp -r ../athena $(BUILD_DIR)/resources/home/dnanexus/

get_python_packages:
	@printf "\n\nDownloading Python dependencies"
	pip download --requirement ../requirements.txt --dest $(BUILD_DIR)/resources/home/dnanexus/packages \
		|| (echo "Failed to download Python packages" && exit 1)

get_bedtools_binary:
	@printf "\n\nDownloading bedtools binary"
	curl -L -o $(BUILD_DIR)/resources/home/dnanexus/bedtools \
		https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
		|| (echo "Failed to download bedtools binary" && exit 1)

build_app:
	@printf "\n\nBuilding DNAnexus app..."
	@dx build --app $(BUILD_DIR)  || (echo "Failed to build DNAnexus app" && exit 1)

clean_up:
	@printf "\n\nCleaning up build dir"
	@rm -r $(BUILD_DIR)
	@printf "\n\nDone.\n"