import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
//Table
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
//Scrollbar
import { Scrollbars } from 'react-custom-scrollbars';
//
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import InputLabel from '@material-ui/core/InputLabel';

const TrackingResultTable = withStyles(theme => ({
	head:{
		backgroundColor: theme.palette.common.black,
		color: theme.palette.common.white,
		position: "sticky",
		top: 0
	},
	body:{
		fontSize: 14,
	},
}))(TableCell);

const styles = theme => ({
	root:{
		width: '100%',
		marginTop: theme.spacing.unit,
	},
	table:{
		minWidth: 10,
	},
	row:{
		'&:nth-of-type(odd)': {
			backgroundColor: theme.palette.background.default,
		},
	},
	formControl: {
	    marginLeft: theme.spacing.unit * 5,
	    marginTop: theme.spacing.unit * 3,
	    minWidth: 240,
	},
	selectEmpty: {
		marginTop: theme.spacing.unit * 2,
	},
});

class Tracking_result extends React.Component {

	constructor(props) {
		super(props);
		this.state = { rownumber:''};
		this.query_track_result = this.query_track_result.bind(this);
	};

	query_track_result(){
		if(this.state.tracking_result == undefined){
			fetch('api/tracking/results/' + window.trackingID, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                tracking_result: result.json,
                tracking_result_shown: result.json })));
		}else{
			clearInterval(this.interval);
		}
	}

	componentDidMount(){
		this.query_track_result(); //delete after testing
		this.interval = setInterval(this.query_track_result, 10000);
	}

	handleChange(event){
		this.setState(state => ({ [event.target.name]: event.target.value }));
		let i = 0;
		let tmp = [];
		for(i; i < event.target.value; i++){
			tmp.push(this.state.tracking_result[i])
		};
		this.setState(state => ({ tracking_result_shown: tmp }));
	};

	turn_on_Tabs(){
		window.tabSwitch = false;
	}

    render() {
    	const { classes } = this.props;
    	const trackResult = this.state.tracking_result_shown;

    	if(this.state.tracking_result == undefined){
    		return(
    			<div>
				<br />
				<br />
				<br />
				<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
					<font> Please hold on ... </font>
				</div>
				<br />
				<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
					<img src={require('./static/waiting.svg')} />
				</div>
				<br />
				<br />
				<br />
			</div>
		);
    	
    	}else{
    		return (
				<div>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Paper className={classes.root}>
							<form autoComplete="off">
								<FormControl className={classes.formControl}>
									<InputLabel>Display rows (Default:100)</InputLabel>
										<Select
										value={this.state.rownumber}
										onChange={this.handleChange.bind(this)}
										name="rownumber"
										className={classes.selectEmpty}
										>
										<MenuItem value={100}>100</MenuItem>
										<MenuItem value={50}>50</MenuItem>
										<MenuItem value={20}>20</MenuItem>
										<MenuItem value={10}>10</MenuItem>
										</Select>
								</FormControl>
							</form>
							<Scrollbars 
							style={{ width: '93%', height: 600, margin:30}}>
								<Table className={classes.table}>
									<TableHead>
										<TableRow>
											<TrackingResultTable align="right">Difference(loci)</TrackingResultTable>
											<TrackingResultTable align="right">BioSample</TrackingResultTable>
											<TrackingResultTable align="right">Strain ID/Alias</TrackingResultTable>
											<TrackingResultTable align="right">SourceSeq</TrackingResultTable>
											<TrackingResultTable align="right">Source country</TrackingResultTable>
											<TrackingResultTable align="right">Isolated year</TrackingResultTable>
											<TrackingResultTable align="right">ST</TrackingResultTable>
											<TrackingResultTable align="right">Serogroup_type</TrackingResultTable>
											<TrackingResultTable align="right">Number of void loci</TrackingResultTable>
										</TableRow>
									</TableHead>
									<TableBody>
										{trackResult.map(row => (
											<TableRow className={classes.row} key={row.BioSample}>
												<TrackingResultTable align="right">{row.distance}</TrackingResultTable>
												<TrackingResultTable align="right">
													<a href={'https://www.ncbi.nlm.nih.gov/biosample/?term='+row.BioSample}
													target="_blank">
														{row.BioSample}
													</a>
												</TrackingResultTable>
												<TrackingResultTable align="right">{row.Strain}</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceSeq}</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceCountry}</TrackingResultTable>
												<TrackingResultTable align="right">{row.IsolatYear}</TrackingResultTable>
												<TrackingResultTable align="right">{row.ST}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Serogroup_type}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Number_of_void_loci}</TrackingResultTable>
											</TableRow>
										))}
									</TableBody>
								</Table>
							</Scrollbars>
							<br />
						</Paper>
					</div>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Link to="/tracking" style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default" onClick={this.turn_on_Tabs}>
								<ReplyIcon />
								&nbsp;&nbsp;
								Back
							</Button>
						</Link>
					</div>
					<br />
				</div>
			);
		}
    }
}

export default withStyles(styles)(Tracking_result);
