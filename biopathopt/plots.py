from typing import Sequence, Optional
import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure

def cdf(
    values: Sequence[float],
    title: str = "Empirical CDF",
    x_label: str = "Value",
    y_label: str = "Cumulative probability",
    height: int = 500,
    width: int = 800,
) -> Figure:
    """Empirical CDF

    This computes the empirical cumulative distribution function:
        F(x) = (number of observations ≤ x) / n

    Args:
        values: A 1D list/array-like of numeric values.
        title: Plot title.
        x_label: X-axis label.
        y_label: Y-axis label.

    Returns:
        A Plotly Figure object with a step CDF curve in the range [0, 1].

    Example:
        >>> data = [3, 1, 2, 2, 5, 4, 4, 4]
        >>> fig = plot_ecdf(data, title="My CDF", rug=True)
        >>> fig.show()
    """
    # Convert to NumPy array of floats and drop NaN/inf values
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("No valid numeric values provided.")

    # Sort values
    x_sorted = np.sort(x)
    n = x_sorted.size

    # Empirical CDF at each observation
    cdf = np.arange(1, n + 1) / n
    # Collapse duplicates so the step plot looks clean:
    # for each unique x, keep the maximum CDF attained at that x
    df = pd.DataFrame({"x": x_sorted, "cdf": cdf})
    df = df.groupby("x", as_index=False)["cdf"].max()

    # Build step plot using a horizontal-then-vertical line shape
    fig = px.line(
        df, 
        x="x", 
        y="cdf", 
        title=title, 
        labels={"x": x_label, "cdf": y_label}, 
        log_x=True,
        width=width,
        height=height
    )
    fig.update_traces(line_shape="hv")  # step function look
    fig.update_yaxes(range=[0, 1], constrain="domain")  # keep CDF within [0,1]
    fig.update_yaxes(range=[0, 1])  # keep CDF within [0,1]
    fig.update_layout(margin=dict(l=40, r=20, t=60, b=40))
    fig.update_xaxes(exponentformat = 'power')
    fig.add_hline(y=0.5, line_width=3, line_dash="dash", line_color="gray")
    return fig

def predicted_v_sim(
    df,
    measured_col,
    sim_col,
    name_col,
    width = 1200,
    height = 800,
    title: str = "Predicted vs Measured",
):
    # Determine plot limits
    min_val = min(df[measured_col].min(), df[sim_col].min())
    max_val = max(df[measured_col].max(), df[sim_col].max())

    # Create scatter plot
    fig = px.scatter(
        df,
        x=measured_col,
        y=sim_col,
        title=title,
        color=name_col,
        text=name_col,
        labels={measured_col: "Measured", sim_col: "Predicted", name_col: "Substrate"},
        width=width,
        height=height,
    )

    # Add perfect correlation line
    fig.add_shape(
        type="line",
        x0=min_val, y0=min_val,
        x1=max_val, y1=max_val,
        line=dict(color="gray", dash="dash"),
        name="Perfect correlation",

    )

    # Adjust text placement and make markers bigger
    fig.update_traces(
        textposition="top center",
        marker=dict(size=14, opacity=0.85)
    )

    #fig.update_layout(template="plotly_white")
    return fig
