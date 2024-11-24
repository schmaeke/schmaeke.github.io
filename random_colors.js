document.addEventListener( "DOMContentLoaded", ( ) => {

  const colors = [
		"#C5233F", "#A11E35",
		"#E2BB36", "#B3962B",
		"#29A245", "#22863A",
		"#58AA9A", "#478B80",
		"#2A8AD9", "#2272B3",
		"#9D3BD6", "#812FB4",
	];

  const hover_colors = new Array( colors.length );

	for ( let i = colors.length - 1; i > 0; i-- )
	{
    const j = Math.floor( Math.random( ) * ( i + 1 ) );
    [ colors[ i ], colors[ j ] ] = [ colors[ j ], colors[ i ] ];
    [ hover_colors[ i ], hover_colors[ j ] ] = [ hover_colors[ j ], hover_colors[ i ] ];
  }

	const amount = 25;
	for ( let i = 0; i < colors.length; i++ )
	{
		hover_colors[ i ] = '#' + colors[ i ].replace( /^#/, '' ).replace( /../g, color => ( '0' + Math.min( 255, Math.max( 0, parseInt( color, 16 ) + amount ) ).toString( 16 ) ).substr( -2 ) );
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
