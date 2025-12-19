// #import "presentation_template.typ": theme
#import "@local/dutten-slides:0.1.0": dutten-theme, progress-marker-both, progress-marker-lower, progress-marker-upper
#import "@local/dutten-common:0.1.0": colors-rgb
#import "@preview/touying:0.6.1": *

#let setup(hw_num, doc) = [
    #show: dutten-theme.with(
    title: [02624 Inverse Problems and Imaging - HW#hw_num Presentation],
    author: [Jeppe Klitgaard],
    date: [2025-12-19],
    institution: [Technical University of Denmark],
    )

    #set block(above: 0.5em, below: 0.5em)

    #doc
]

#let color-func(color, c) = [
    #set text(fill: color)
    #c
]
#let blue = color-func.with(colors-rgb.at("blue"))
#let purple = color-func.with(colors-rgb.at("purple"))
#let red = color-func.with(colors-rgb.at("corporate red"))
#let side-by-side(columns: (1fr, 1fr), left, right) = grid(
columns: columns,
column-gutter: 1em,
left, right
)
