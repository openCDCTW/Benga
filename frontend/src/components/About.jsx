import React from 'react';
import ReactDOM from 'react-dom';
import Paper from '@material-ui/core/Paper';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';

const styles = theme => ({
	paper_root: {
		...theme.mixins.gutters(),
		paddingTop: theme.spacing.unit * 2,
		paddingBottom: theme.spacing.unit * 3,
		width:'90%',
	},
})

class About extends React.Component {
	render() {
		const { classes } = this.props;

    	return (
    		<div>
    			<br />
    			<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
	    			<Paper className={classes.paper_root} elevation={4}>
	    				<Typography component="p" style={{ width:'90%',textAlign: 'justify',
							margin:'auto',fontSize:16 }}>
							
						</Typography>
	    			</Paper>
    			</div>
    			<br />
    		</div>
    	)
    }
}

export default withStyles(styles)(About)