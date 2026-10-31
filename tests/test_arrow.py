import pytest
import numpy as np

from homologyviz.arrow import Arrow


@pytest.mark.parametrize(
    "ratio",
    [0.01, 0.1, 0.5, 0.99, 1.0],
)
def test_valid_ratio_tail_head_width(ratio):
    arrow = Arrow(
        x1=100,
        x2=1000,
        y=10,
        ratio_tail_head_width=ratio,
    )

    assert arrow.ratio_tail_head_width == ratio
    assert arrow.tail_width == pytest.approx(arrow.head_width * ratio)


@pytest.mark.parametrize(
    "ratio",
    [-1.0, -0.5, 0.0, 1.01, 1.5, 2.0],
)
def test_invalid_ratio_tail_head_width(ratio):
    with pytest.raises(ValueError, match="ratio_tail_head_width"):
        Arrow(
            x1=100,
            x2=1000,
            y=10,
            ratio_tail_head_width=ratio,
        )


@pytest.mark.parametrize(
    "head_width",
    [-2, -1, 0],
)
def test_invalid_head_width(head_width):
    with pytest.raises(ValueError, match="head_width"):
        Arrow(
            x1=100,
            x2=1000,
            y=10,
            head_width=head_width,
        )


@pytest.mark.parametrize(
    "head_height",
    [-200, -1, 0],
)
def test_invalid_head_height(head_height):
    with pytest.raises(ValueError, match="head_height"):
        Arrow(
            x1=100,
            x2=1000,
            y=10,
            head_height=head_height,
        )


def test_identical_coordinates():
    with pytest.raises(ValueError, match="x1 and x2"):
        Arrow(
            x1=500,
            x2=500,
            y=10,
        )


def test_arrow_repr():
    arrow = Arrow(
        x1=100,
        x2=1000,
        y=10,
        ratio_tail_head_width=0.5,
        head_width=2,
        head_height=200,
    )

    expected = (
        "Arrow(x1=100, x2=1000, y=10, "
        "ratio_tail_head_width=0.5, "
        "head_width=2, head_height=200)"
    )

    assert repr(arrow) == expected


def test_forward_arrow_tip():
    arrow = Arrow(x1=100, x2=1000, y=10)

    x_values, y_values = arrow.get_coordinates()

    assert x_values[4] == 1000
    assert y_values[4] == 10


def test_forward_arrow_coordinates():
    arrow = Arrow(
        x1=100,
        x2=1000,
        y=10,
        ratio_tail_head_width=0.5,
        head_width=2,
        head_height=200,
    )

    x_values, y_values = arrow.get_coordinates()

    expected_x = np.array([100, 100, 800, 800, 1000, 800, 800, 100, 100])

    expected_y = np.array([10, 10.5, 10.5, 11, 10, 9, 9.5, 9.5, 10])

    np.testing.assert_allclose(x_values, expected_x)
    np.testing.assert_allclose(y_values, expected_y)


def test_reverse_arrow_coordinates():
    arrow = Arrow(
        x1=1000,
        x2=100,
        y=10,
        ratio_tail_head_width=0.5,
        head_width=2,
        head_height=200,
    )

    x_values, y_values = arrow.get_coordinates()

    expected_x = np.array([1000, 1000, 300, 300, 100, 300, 300, 1000, 1000])

    expected_y = np.array([10, 9.5, 9.5, 9, 10, 11, 10.5, 10.5, 10])

    np.testing.assert_allclose(x_values, expected_x)
    np.testing.assert_allclose(y_values, expected_y)


@pytest.mark.parametrize(
    "x1, x2",
    [
        (100, 150),  # Short forward arrow
        (150, 100),  # Short reverse arrow
        (100, 300),  # Forward arrow exactly equal to head_height
        (300, 100),  # Reverse arrow exactly equal to head_height
    ],
)
def test_short_arrow_coordinates(x1, x2):
    arrow = Arrow(
        x1=x1,
        x2=x2,
        y=10,
        head_width=2,
        head_height=200,
    )

    x_values, y_values = arrow.get_coordinates()

    assert len(x_values) == 5
    assert len(y_values) == 5

    # The arrow tip must always be at (x2, y).
    assert x_values[2] == x2
    assert y_values[2] == 10

    # The polygon must start and end at the same point.
    assert x_values[0] == x_values[-1]
    assert y_values[0] == y_values[-1]


def test_forward_reverse_symmetry():
    x1 = 100
    x2 = 1000
    y = 10

    forward = Arrow(x1=x1, x2=x2, y=y)
    reverse = Arrow(x1=x2, x2=x1, y=y)

    forward_x, forward_y = forward.get_coordinates()
    reverse_x, reverse_y = reverse.get_coordinates()

    expected_reverse_x = x1 + x2 - forward_x
    expected_reverse_y = 2 * y - forward_y

    np.testing.assert_allclose(reverse_x, expected_reverse_x)
    np.testing.assert_allclose(reverse_y, expected_reverse_y)


def test_arrow_properties():
    arrow = Arrow(
        x1=100,
        x2=1000,
        y=10,
        ratio_tail_head_width=0.5,
        head_width=2,
    )

    assert arrow.tail_width == 1

    # Modify head width.
    arrow.head_width = 4

    assert arrow.tail_width == 2

    # Modify the ratio.
    arrow.ratio_tail_head_width = 0.25

    assert arrow.tail_width == 1


if __name__ == "__main__":
    test_forward_reverse_symmetry()
