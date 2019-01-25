import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import { Link } from 'react-router-dom';
import { withStyles } from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';

const styles = theme => ({
    textField:{
        margin: theme.spacing.unit * 3,
        width:'95%'
    },
    container:{
        display:'flex',
        flexWrap:'wrap',
    },
})

class TrackingSearch extends React.Component {
	constructor(props) {
		super(props);
		this.state = {};
	}

    _onKeyPress(event){
	    if(event.charCode === 13){
	        event.preventDefault();
	    }
    }

    BioSample_handleChange(event){
    	this.setState(state => ({ BioSampleID: this.search.value }));
    }

	render(){

		const { classes } = this.props;

		return(
			<div>
				<br />
				<br />
	            <Paper>
	                <div>
	                    <form className={classes.container}>
	                        <TextField
	                            inputRef={ID => this.search = ID}
	                            label="BioSample ID"
	                            type="search"
	                            placeholder="Input BioSample ID here"
	                            onChange={this.BioSample_handleChange.bind(this)}
	                            className={classes.textField}
	                            onKeyPress={this._onKeyPress}
	                            margin="normal"
	                            variant="outlined"
	                            helperText=""
	                        />
	                    </form>
	                </div>
	                <br />
	                <div style={{ width:'97%', display:'flex', justifyContent:'flex-end', 
	                        alignItems:'flex-end'}}>
	                    <Button variant="contained" color="default">
	                        Search
	                    </Button>
	                </div>
	                <br />
	            </Paper>
	            <br />
	            <br />
            </div>
    	);
	}
}

export default withStyles(styles)(TrackingSearch);