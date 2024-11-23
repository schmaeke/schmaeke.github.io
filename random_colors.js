document.addEventListener( "DOMContentLoaded", ( ) => {

  const colors = [
		"#C5233F", "#DC2C48",
		"#E2BB36", "#E9C946",
		"#29A245", "#33B34E",
		"#58AA9A", "#66B8A7",
		"#2A8AD9", "#3C9BE1",
		"#9D3BD6", "#A948E0",
	];

  const hover_colors = [
		"#DC2C48", "#E94C66",
		"#E9C946", "#F3D870",
		"#33B34E", "#42CC60",
		"#66B8A7", "#80C9BC",
		"#3C9BE1", "#5FB2F0",
		"#A948E0", "#BB6EF1",
	];

	for ( let i = colors.length - 1; i > 0; i-- )
	{
    const j = Math.floor( Math.random( ) * ( i + 1 ) );
    [ colors[ i ], colors[ j ] ] = [ colors[ j ], colors[ i ] ];
    [ hover_colors[ i ], hover_colors[ j ] ] = [ hover_colors[ j ], hover_colors[ i ] ];
  }

  const links = document.querySelectorAll( "main > a" );

  links.forEach( ( link, index ) => {
		const color_index = index % colors.length;
		const base_color = colors[ color_index ];
		const hover_color = hover_colors[ color_index ];

    link.style.backgroundColor = base_color;

		// Set hover effect
		link.addEventListener("mouseover", () => {
			link.style.backgroundColor = hover_color;
		});

		link.addEventListener("mouseout", () => {
			link.style.backgroundColor = base_color;
		});
  } );
} );
