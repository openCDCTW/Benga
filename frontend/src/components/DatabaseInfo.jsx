import React from 'react';
import ReactDOM from 'react-dom';
import Typography from '@material-ui/core/Typography';

export default class DatabaseInfo extends React.Component {

	constructor(props) {
		super(props);
	}

	render(){

		const { classes } = this.props;

		if(this.props.database == 1){
			return (
				<div>
					<br />
					<Typography component="p" style={{ width:'90%',textAlign: 'justify',
					margin:'auto',fontSize:16, marginTop: '60px' }}>
					Core gene: 2951
					</Typography>
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	                	<img style={{ width:'95%'}} 
	                	src={require('./static/Vibrio_cholerae_locus_distribution.png')} />
	            	</div>
				</div>
			)
		}
	}

}