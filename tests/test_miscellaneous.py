import pytest
import tempfile
from pathlib import Path
import platform
from unittest.mock import patch, MagicMock

from homologyviz.logger import get_logger
from homologyviz.miscellaneous import (
    get_package_path,
    get_os,
    is_blastn_installed,
    round_up_to_nearest_significant_digit,
    delete_files,
    clean_directory,
)


# ==== Test get_package_path()
def test_get_package_path():
    package_path = get_package_path("homologyviz")
    assert isinstance(package_path, Path)
    assert package_path.exists()
    assert package_path.is_dir()


# ==== Test get_os()
def test_get_os_windows():
    with patch.object(platform, "system", return_value="Windows"):
        assert get_os() == "Windows"


def test_get_os_macos():
    with patch.object(platform, "system", return_value="Darwin"):
        assert get_os() == "macOS"


def test_get_os_linux():
    with patch.object(platform, "system", return_value="Linux"):
        assert get_os() == "Linux"


def test_get_os_unknown():
    with patch.object(platform, "system", return_value="HaikuOS"):
        assert get_os() == "Unknown"


# ==== Test is_blatn_installed()
def test_is_blatn_installed_success():
    with patch("subprocess.run") as mock_run:
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "blastn: 2.12.0"
        mock_run.return_value = mock_result

        assert is_blastn_installed() is True


def test_is_blatn_installed_not_found():
    with patch("subprocess.run", side_effect=FileNotFoundError):
        assert is_blastn_installed() is False


def test_is_blatn_installed_exception():
    with patch("subprocess.run", side_effect=Exception):
        assert is_blastn_installed() is False


def test_is_blastn_installed_error():
    with patch("subprocess.run") as mock_run:
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "some error"
        mock_run.return_value = mock_result

        assert is_blastn_installed() is False


# ==== Test delete_files()
# def test_delete_files_logs(tmp_path, caplog):
#     # Create a real file that will be deleted
#     existing_file = tmp_path / "file1.txt"
#     existing_file.write_text("hello")

#     # Create a fake file that doesn't exist
#     missing_file = tmp_path / "file2.txt"

#     logger = get_logger("homologyviz.miscellaneous")

#     # Capture logs at INFO level and higher
#     with caplog.at_level("INFO", logger="homologyviz.miscellaneous"):
#         delete_files([existing_file, missing_file])

#     # Assert log messages
#     logs = caplog.text
#     assert f"Deleted file: {existing_file}" in logs
#     assert f"File not found: {missing_file}" in logs
#     assert not existing_file.exists()  # Ensure it was actually deleted

# def test_delete_files_logs(tmp_path, caplog):
#     # Create test files
#     existing_file = tmp_path / "file1.txt"
#     existing_file.write_text("hello")
#     missing_file = tmp_path / "file2.txt"

#     # Get the logger used in the module being tested
#     logger = get_logger("homologyviz.miscellaneous")

#     # Save original propagate value
#     original_propagate = logger.propagate
#     logger.propagate = (
#         True  # Allow logs to bubble up to root so pytest can capture them
#     )

#     try:
#         with caplog.at_level("INFO"):
#             delete_files([existing_file, missing_file])
#     finally:
#         logger.propagate = original_propagate  # Restore to avoid affecting other tests

#     logs = caplog.text
#     assert f"Deleted file: {existing_file}" in logs
#     assert f"File not found: {missing_file}" in logs
#     assert not existing_file.exists()


# TODO fix this test, it is not working
def test_delete_files_logs(tmp_path, caplog):
    # Create test files
    existing_file = tmp_path / "file1.txt"
    existing_file.write_text("hello")
    missing_file = tmp_path / "file2.txt"

    # Get the logger used in delete_files
    logger = get_logger("homologyviz.miscellaneous")

    # Backup and remove custom handlers temporarily
    original_handlers = logger.handlers[:]
    for handler in original_handlers:
        logger.removeHandler(handler)

    # Enable propagation so logs bubble to the root logger
    original_propagate = logger.propagate
    logger.propagate = True

    try:
        with caplog.at_level("INFO"):
            delete_files([existing_file, missing_file])
    finally:
        # Restore logger state
        logger.handlers = original_handlers
        logger.propagate = original_propagate

    logs = caplog.text
    assert f"Deleted file: {existing_file}" in logs
    assert f"File not found: {missing_file}" in logs
    assert not existing_file.exists()


# ==== Test clean_directory()
def test_clean_directory():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        file1 = Path(temp_dir) / "test1.txt"
        file2 = Path(temp_dir) / "test2.txt"
        file1.write_text("lol")
        file2.write_text("dummy")
        clean_directory(temp_dir)
        assert list(temp_dir.iterdir()) == []


def test_clean_directory_is_empty():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        assert clean_directory(temp_dir) is None


# ==== Test round_up_to_nearest_significant_digit()
def test_round_up_to_nearest_significant_digit():
    assert round_up_to_nearest_significant_digit(153) == 200
    assert round_up_to_nearest_significant_digit(89) == 90
    assert round_up_to_nearest_significant_digit(1000) == 1000
    assert round_up_to_nearest_significant_digit(1) == 1


def test_round_up_to_nearest_significant_digit_invalid():
    with pytest.raises(ValueError):
        round_up_to_nearest_significant_digit(0)
    with pytest.raises(ValueError):
        round_up_to_nearest_significant_digit(-10)
