from pathlib import Path
import base64
import binascii


def save_uploaded_file(
    file_name: str,
    content: str,
    temp_folder_path: Path,
) -> str | None:
    """Decode the content and write it to a temporary file.

    Returns the file path as a string if successful, otherwise returns None.
    """
    try:
        # Ensure content is properly formatted
        if ";base64," not in content:
            raise ValueError("Content is not base64-encoded or improperly formatted.")

        # Decode content
        data = content.split(";base64,")[1]
        decoded_data = base64.b64decode(data)

        # Save uploaded file
        output_path = temp_folder_path / file_name
        with open(output_path, "wb") as f:
            f.write(decoded_data)

        # Dash doesn't like Path; hence, we need to cast Path to str.
        return str(output_path)

    except (IndexError, ValueError, binascii.Error) as e:
        print(f"Failed to decode and save uplaoded file: {e}")
        return None
