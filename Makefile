.PHONY: all clean

# ビルドディレクトリの設定
BUILD_DIR := .build

# デフォルトターゲット
all:
	@echo "=== Configuring and Building in $(BUILD_DIR) ==="
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(shell pwd) ..
	@cd $(BUILD_DIR) && make -j$(shell nproc)
	@cd $(BUILD_DIR) && make install
	@echo "=== Build Complete ==="
	@echo "Executable is located in bin/KVCOpticalSim"

# クリーンターゲット（ビルドディレクトリとbinを削除）
clean:
	@echo "=== Cleaning build artifacts ==="
	@rm -rf $(BUILD_DIR) bin
	@echo "=== Clean Complete ==="
