import pytest

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
