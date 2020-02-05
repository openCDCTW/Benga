import React from 'react';
import ReactDOM from 'react-dom';
import Paper from '@material-ui/core/Paper';
import { withStyles } from '@material-ui/core/styles';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
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
    selectEmpty: {
    	marginTop: theme.spacing.unit * 2,
	},
	formControl: {
	    margin: theme.spacing.unit,
	    minWidth: 200,
	},
})

class About extends React.Component {

	constructor(props) {
		super(props);
        this.state = { databaseInfo: 'VC', };
	}

	databaseInfo_Change(event){
		this.setState(state => ({
			[event.target.name]: event.target.value
		}));
	}

	render() {
		const { classes } = this.props;

    	return (
    		<div>
    			<br />
    			<div style={{ justifyContent:'center',margin:'auto',display:'flex'}}>
	    			<Paper className={classes.paper_root} elevation={4}>
	    				<h2 style={{ margin:'30px' }}>Database Information</h2>
	    				<div style={{ float:'left', marginLeft:'25px' }}>
	                        <form className={classes.selectcss} autoComplete="off">
	                            <FormControl required className={classes.formControl} >
	                              <InputLabel htmlFor="database-required">Database</InputLabel>
	                                <Select
	                                  value={this.state.databaseInfo}
	                                  onChange={this.databaseInfo_Change.bind(this)}
	                                  name="databaseInfo"
	                                  inputProps={{
	                                    id: 'database-required',
	                                  }}
	                                  className={classes.selectEmpty}
	                                  >
	                                  <MenuItem value={'VC'}>Vibrio cholerae</MenuItem>
	                                  <MenuItem value={'NM'}>Neisseria meningitidis</MenuItem>
	                                </Select>
	                            </FormControl>
	                        </form>
                    	</div>
                    	<div style={{ marginTop:'110px'}}>
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