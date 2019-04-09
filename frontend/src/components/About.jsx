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

	handleChange(event){
        this.setState(state => ({ [event.target.name]: event.target.value }));
    }

	render() {
		const { classes } = this.props;

    	return (
    		<div>
    			<br />
    			<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
	    			<Paper className={classes.paper_root} elevation={4}>
	    				<h2>&nbsp;Database Information</h2>
	    				<div style={{ float:'left', marginLeft:'15px' }}>
	                        <form className={classes.selectcss} autoComplete="off">
	                            <FormControl required className={classes.formControl} >
	                              <InputLabel>Species</InputLabel>
	                                <Select
	                                  value={this.state.databaseInfo}
	                                  onChange={this.handleChange.bind(this)}
	                                  name="databaseInfo"
	                                  >
	                                  <MenuItem value={1}>Vibrio cholerae</MenuItem>
	                                </Select>
	                            </FormControl>
	                        </form>
                        </div>
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