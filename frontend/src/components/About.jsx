import React from 'react';
import ReactDOM from 'react-dom';
import Paper from '@material-ui/core/Paper';
import { withStyles } from '@material-ui/core/styles';
// import Typography from '@material-ui/core/Typography';
//
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
//
import DatabaseInfo from './DatabaseInfo.jsx';



const styles = theme => ({
	paper_root: {
		...theme.mixins.gutters(),
		paddingTop: theme.spacing.unit * 2,
		paddingBottom: theme.spacing.unit * 3,
		width:'90%',
	},
	selectcss: {
	    display: 'flex',
	    flexWrap: 'wrap',
	},
	formControl: {
	    margin: theme.spacing.unit,
	    minWidth: 200,
	},
})

class About extends React.Component {

	constructor(props) {
		super(props);
		this.state = {
			databaseInfo: 1,
		};
	}


	render() {
		const { classes } = this.props;

    	return (
    		<div>
    			<br />
    			<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
	    			<Paper className={classes.paper_root} elevation={4}>
	    				<h2 style={{ margin:'30px' }}>Database Information</h2>
                    	<br />
                    	<div>
                    		<DatabaseInfo database={this.state.databaseInfo}/>
                    	</div>
	    			</Paper>
    			</div>
    			<br />
    		</div>
    	)
    }
}

export default withStyles(styles)(About)