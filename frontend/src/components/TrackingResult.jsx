import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
//Table
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';

const TrackingResultTable = withStyles(theme => ({
	head:{
		backgroundColor: theme.palette.common.black,
		color: theme.palette.common.white,
	},
	body:{
		fontSize: 14,
	},
}))(TableCell);

const styles = theme => ({
	root:{
		width: '100%',
		margibTop: theme.spacing.unit * 5,
		overflowX: 'auto',
	},
	table:{
		minWidth: 10,
	},
	row:{
		'&:nth-of-type(odd)': {
			backgroundColor: theme.palette.background.default,
		},
	}
});

class Tracking_result extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_track_result = this.query_track_result.bind(this);
	};

	query_track_result(){
		if(this.state.tracking_result == undefined){
			fetch('api/tracking/results/' + window.trackingID, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                tracking_result: result.json })));
		}else{
			clearInterval(this.interval);
		}
	}
	componentDidMount(){
		this.interval = setInterval(this.query_track_result, 10000);
	}

    render() {
    	const { classes } = this.props;
    	const trackResult = this.state.tracking_result;

    	if(this.state.tracking_result == undefined){
    		return(
    			<div>
				<Paper square>
					<Tabs value={false} centered>
						<Tab label=" " disabled/>
					</Tabs>
				</Paper>
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
				<div id="url">
					<Paper square>
						<Tabs value={false} centered>
							<Tab label=" " disabled/>
						</Tabs>
					</Paper>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Paper className={classes.root}>
							<Table className={classes.table}>
								<TableHead>
									<TableRow>
										<TrackingResultTable align="right">Distance</TrackingResultTable>
										<TrackingResultTable align="right">BioSample</TrackingResultTable>
										<TrackingResultTable align="right">Country</TrackingResultTable>
										<TrackingResultTable align="right">ST</TrackingResultTable>
										<TrackingResultTable align="right">Serogroup_serotype</TrackingResultTable>
										<TrackingResultTable align="right">Strain</TrackingResultTable>
										<TrackingResultTable align="right">Year</TrackingResultTable>
									</TableRow>
								</TableHead>
								<TableBody>
									{trackResult.map(row => (
										<TableRow className={classes.row} key={row.BioSample}>
											<TrackingResultTable align="right">{row.distance}</TrackingResultTable>
											<TrackingResultTable align="right">{row.BioSample}</TrackingResultTable>
											<TrackingResultTable align="right">{row.Country}</TrackingResultTable>
											<TrackingResultTable align="right">{row.ST}</TrackingResultTable>
											<TrackingResultTable align="right">{row.Serogroup_serotype}</TrackingResultTable>
											<TrackingResultTable align="right">{row.Strain}</TrackingResultTable>
											<TrackingResultTable align="right">{row.Year}</TrackingResultTable>
										</TableRow>
									))}
								</TableBody>
							</Table>
						</Paper>
					</div>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Link to="/tracking" style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">
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