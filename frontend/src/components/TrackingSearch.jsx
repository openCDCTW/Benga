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
import SearchIcon from '@material-ui/icons/Search';
import Typography from '@material-ui/core/Typography';

import OutlinedInput from '@material-ui/core/OutlinedInput';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';


const styles = theme => ({
    textField:{
        marginLeft: '23px',
        marginTop: '25px',
        width:'30%'
    },
    container:{
        display:'flex',
        flexWrap:'wrap',
    },
    numberField:{
    	marginLeft: '10px',
        marginTop: '10px',
        marginRight: '10px',
        width:'110px'
    },
    countrySelect:{
    	marginLeft: '23px',
        marginTop: '25px',
    	width: '30%',
    }
})

class TrackingSearch extends React.Component {
	constructor(props) {
		super(props);

		let nowYear = new Date();

		this.state = { 
			yearError: false,
			country: '',
			labelWidth: 0,
			nowYear: nowYear.getFullYear(),
		};
	}

	componentDidMount() {
	    this.setState({
	      labelWidth: ReactDOM.findDOMNode(this.InputLabelRef).offsetWidth,
	    });
	}

    _onKeyPress(event){
	    if(event.charCode === 13){
	        event.preventDefault();
	    }
    }

    bioSampleHandleChange(event){
    	this.setState(state => ({ bioSampleID: this.biosample.value }));
    }

    strainHandleChange(event){
    	this.setState(state => ({ strain: this.strain.value }));
    }

    yearFromHandleChange(event){
    	let year_from = Number.parseInt(this.yearFrom.value, 10);

    	if( this.yearFrom.value == year_from ){
    		if( year_from > 0 && year_from <= this.state.nowYear ){
	    		this.setState(state => ({ 
	    			yearFrom: this.yearFrom.value , 
	    			yearError:false
	    		}));
	    	}else{
	    		this.yearFrom.value = '';
	    		this.setState(state => ({ yearFrom: undefined }));
	    	}
    	}else{
    		this.yearFrom.value = '';
    		this.setState(state => ({ yearFrom: undefined }));
    	}
    }

    yearToHandleChange(event){
    	let year_To = Number.parseInt(this.yearTo.value, 10);
    	
    	if( this.yearTo.value == year_To ){
    		if( year_To > 0 && year_To <= this.state.nowYear ){
	    		this.setState(state => ({ 
	    			yearTo: this.yearTo.value, 
	    			yearError: false
	    		}));
	    	}else{
	    		this.yearTo.value = '';
	    		this.setState(state => ({ yearTo: undefined }));
	    	}
    	}else{
    		this.yearTo.value = '';
    		this.setState(state => ({ yearTo: undefined }));
    	}
    }

    countryHandleChange(event){
    	this.setState(state => ({ [event.target.name]: event.target.value }));
    }

    search(){
    	console.log(this.state.bioSampleID);
    	console.log(this.state.strain);
    	console.log(Number.parseInt(this.state.yearFrom, 10));
    	console.log(Number.parseInt(this.state.yearTo, 10));

    	if( this.state.yearFrom > this.state.yearTo ){
    		this.setState(state => ({ yearError: true }));
    		alert('Input year error');
    	}
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
	                            inputRef={ID => this.biosample = ID}
	                            label="BioSample ID"
	                            type="search"
	                            placeholder="Input BioSample ID here"
	                            onChange={this.bioSampleHandleChange.bind(this)}
	                            className={classes.textField}
	                            onKeyPress={this._onKeyPress}
	                            margin="normal"
	                            variant="outlined"
	                        />
	                        <TextField
	                            inputRef={ID => this.strain = ID}
	                            label="Strain"
	                            type="search"
	                            placeholder="Input strain here"
	                            onChange={this.strainHandleChange.bind(this)}
	                            className={classes.textField}
	                            onKeyPress={this._onKeyPress}
	                            margin="normal"
	                            variant="outlined"
	                        />
	                        <FormControl variant="outlined" className={classes.countrySelect}>
	                        	<InputLabel
									ref={ref => {
									  this.InputLabelRef = ref;
									}}
								>
								Country
								</InputLabel>
								<Select
								value={this.state.country}
								onChange={this.countryHandleChange.bind(this)}
								name="country"
								input={
										<OutlinedInput
										labelWidth={this.state.labelWidth}
										name="country"
										/>
									}
								MenuProps={{
									PaperProps: {
										style: {
											maxHeight: 49 * 4.5 + 8,
									    },
									  },
								  }}
								>
									<MenuItem value=''>Any</MenuItem>
									<MenuItem value='Taiwan'>Taiwan</MenuItem>
									<MenuItem value='China'>China</MenuItem>
									<MenuItem value='Russia'>Russia</MenuItem>
									<MenuItem value='Ukraine'>Ukraine</MenuItem>
									<MenuItem value='USA'>USA</MenuItem>
								</Select>
							</FormControl>
	                    </form>
	                </div>
	                <div>
	                	<form className={classes.container}>
	                		<Typography style={{ fontSize:'18px', marginLeft:'25px', 
	                		display:'flex', justifyContent:'center', alignItems:'center'}}> 
	                		from </Typography>
	                        &nbsp;&nbsp;
	                        <TextField
	                        	label="Year"
	                			inputRef={year => this.yearFrom = year}
	                            type="number"
	                            onChange={this.yearFromHandleChange.bind(this)}
	                            className={classes.numberField}
	                            inputProps={{ 
	                            	style:{ textAlign: 'center'},
	                            }}
	                            margin="normal"
	                            variant="outlined"
	                            error={this.state.yearError}
	                        />
	                        &nbsp;&nbsp;
	                        <Typography style={{ fontSize:'18px', display:'flex', justifyContent:'center', 
	                        alignItems:'center'}}>to</Typography>
	                        &nbsp;&nbsp;
	                        <TextField
	                        	label="Year"
	                        	inputRef={year => this.yearTo = year}
	                            type="number"
	                            onChange={this.yearToHandleChange.bind(this)}
	                            className={classes.numberField}
	                            inputProps={{ 
	                            	style:{ textAlign: 'center'},
	                            }}
	                            margin="normal"
	                            variant="outlined"
	                            error={this.state.yearError}
	                        />
	                    </form>
	                </div>
	                <br />
	                <div style={{ width:'97%', display:'flex', justifyContent:'flex-end', 
	                        alignItems:'flex-end'}}>
	                    <Button variant="contained" color="default" onClick={this.search.bind(this)}>
	                        Search
	                        &nbsp;&nbsp;
	                        <SearchIcon />
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