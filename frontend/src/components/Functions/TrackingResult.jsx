import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import { Prompt } from 'react-router';
import { withStyles } from '@material-ui/core/styles';
import { Scrollbars } from 'react-custom-scrollbars';
import Paper from '@material-ui/core/Paper';
import Button from '@material-ui/core/Button';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import CircularProgress from '@material-ui/core/CircularProgress';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import InputLabel from '@material-ui/core/InputLabel';
import ReplyIcon from '@material-ui/icons/Reply';
import CloudDownloadIcon from '@material-ui/icons/CloudDownload';
import green from '@material-ui/core/colors/green';

const TrackingResultTable = withStyles(theme => ({
	head:{
		backgroundColor: theme.palette.common.black,
		color: theme.palette.common.white,
		position: "sticky",
		top: 0,
	},
	body:{
		fontSize: 14,
	},
}))(TableCell);

const styles = theme => ({
	root:{
		width: '100%',
		marginTop: theme.spacing(1),
	},
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
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
	    marginLeft: theme.spacing(5),
	    marginTop: theme.spacing(3),
	    minWidth: 240,
	},
	selectEmpty: {
		marginTop: theme.spacing(2),
	},
	downloadButton:{
        color: theme.palette.getContrastText(green[900]),
        backgroundColor: green[600],
        '&:hover': {
            backgroundColor:green[900],
        },
    },
});

class TrackingResult extends React.Component {

	constructor(props) {
		super(props);
		this.state = { rowNumber:''};
		this.queryTrackResult = this.queryTrackResult.bind(this);
	};

	queryTrackResult(){
		if(this.state.trackingResult == undefined){
			fetch('/cgMLST/api/tracking/results/' + window.trackingID + '/', { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                trackingResult: result.json,
                trackingResultShow: result.json,
                trackingZip: result.zip })));
		}else{
			clearInterval(this.interval);
		}
	}

	componentDidMount(){
		this.interval = setInterval(this.queryTrackResult, 3000);
	}

	handleChange(event){
		this.setState(state => ({ [event.target.name]: event.target.value }));
		let i = 0;
		let tmp = [];

		for(i; i < event.target.value; i++){
			tmp.push(this.state.trackingResult[i])
		};

		this.setState(state => ({ trackingResultShow: tmp }));
	};

    render() {
    	const { classes } = this.props;
    	const trackResult = this.state.trackingResultShow;

    	if(this.state.trackingResult == undefined){
    		return(
    			<div style={{ marginTop:'100px' }}>
    				<Prompt
                            when={true}
                            message="Are you sure to leave now?"/>
					<div className={classes.displayCenter}>
						<font size="6"> Please hold on ... </font>
					</div>
					<div className={classes.displayCenter}>
					    <CircularProgress style={{ margin:'50px' }} size={175} />
	                </div>
				</div>
		    );
    	}else if(this.state.trackingResult != undefined && window.database == 'Vibrio_cholerae'){
    		return (
				<div>
					<Prompt
                            when={true}
                            message="You are leaving the page. Please save result, or it will lose. Are you sure to leave now?"/>
					<div className={classes.displayCenter}>
						<Paper className={classes.root}>
							<div>
                                <form autoComplete="off">
                                    <FormControl className={classes.formControl}>
                                        <InputLabel>Display rows (Default:100)</InputLabel>
                                            <Select
                                            value={this.state.rowNumber}
                                            onChange={this.handleChange.bind(this)}
                                            name="rowNumber"
                                            className={classes.selectEmpty}
                                            >
                                            <MenuItem value={100}>100</MenuItem>
                                            <MenuItem value={50}>50</MenuItem>
                                            <MenuItem value={20}>20</MenuItem>
                                            <MenuItem value={10}>10</MenuItem>
                                            </Select>
                                    </FormControl>
                                </form>
								<div style={{ float:'right', marginTop:'-35px', marginRight:'35px' }}>
	                        		<a download href={this.state.trackingZip} style={{ textDecoration:'none' }}>
		                        		<Button variant="contained" className={classes.downloadButton}>
		                                	Download
		                                	&nbsp;&nbsp;
		                                	<CloudDownloadIcon />
		                        		</Button>
	                        		</a>
	                    		</div>
	                    	</div>
							<Scrollbars
							style={{ width: '93%', height: 600, margin:30}}>
								<Table className={classes.table}>
									<TableHead>
										<TableRow>
											<TrackingResultTable align="right">Difference(loci)</TrackingResultTable>
											<TrackingResultTable align="right">BioSample no.</TrackingResultTable>
											<TrackingResultTable align="right">Identifier</TrackingResultTable>
											<TrackingResultTable align="right">Source country</TrackingResultTable>
											<TrackingResultTable align="right">Isolated year</TrackingResultTable>
											<TrackingResultTable align="right">Serogroup type</TrackingResultTable>
											<TrackingResultTable align="right">ST</TrackingResultTable>
											<TrackingResultTable align="right">SourceSeq.</TrackingResultTable>
											<TrackingResultTable align="right">No. void loci</TrackingResultTable>
										</TableRow>
									</TableHead>
									<TableBody>
										{trackResult.map(row => (
											<TableRow className={classes.row} key={row.BioSample}>
												<TrackingResultTable align="right">{row["Distance (loci)"]}</TrackingResultTable>
												<TrackingResultTable align="right">
													<a href={'https://www.ncbi.nlm.nih.gov/biosample/?term='+row.BioSample}
													target="_blank">
														{row.BioSample}
													</a>
												</TrackingResultTable>
												<TrackingResultTable align="right">{row.Identifier}</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceCountry}</TrackingResultTable>
												<TrackingResultTable align="right">{row.IsolatYear}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Serogroup_type}</TrackingResultTable>
												<TrackingResultTable align="right">{row.ST}</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceSeq}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Void_loci}</TrackingResultTable>
											</TableRow>
										))}
									</TableBody>
								</Table>
							</Scrollbars>
						</Paper>
					</div>
					<div className={classes.displayCenter}>
						<Link to="/cgMLST/non-release/tracking" style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">
								<ReplyIcon />
								&nbsp;&nbsp;
								Back
							</Button>
						</Link>
					</div>
				</div>
			);
		}else if(this.state.trackingResult != undefined && window.database == 'Neisseria_meningitidis'){
    		return (
				<div>
					<Prompt
                            when={true}
                            message="You are leaving the page. Please save result, or it will lose. Are you sure to leave now?"/>
					<div className={classes.displayCenter}>
						<Paper className={classes.root}>
							<div>
								<div>
									<form autoComplete="off">
										<FormControl className={classes.formControl}>
											<InputLabel>Display rows (Default:100)</InputLabel>
												<Select
												value={this.state.rowNumber}
												onChange={this.handleChange.bind(this)}
												name="rowNumber"
												className={classes.selectEmpty}
												>
												<MenuItem value={100}>100</MenuItem>
												<MenuItem value={50}>50</MenuItem>
												<MenuItem value={20}>20</MenuItem>
												<MenuItem value={10}>10</MenuItem>
												</Select>
										</FormControl>
									</form>
								</div>
								<div style={{ float:'right', marginTop:'-35px', marginRight:'35px' }}>
	                        		<a download href={this.state.trackingZip} style={{ textDecoration:'none' }}>
		                        		<Button variant="contained" className={classes.downloadButton}>
		                                	Download
		                                	&nbsp;&nbsp;
		                                	<CloudDownloadIcon />
		                        		</Button>
	                        		</a>
	                    		</div>
	                    	</div>
							<Scrollbars
							style={{ width: '93%', height: 600, margin:30}}>
								<Table className={classes.table}>
									<TableHead>
										<TableRow>
                                            <TrackingResultTable align="right">Difference(loci)</TrackingResultTable>
                                            <TrackingResultTable align="right">BioSample no.</TrackingResultTable>
                                            <TrackingResultTable align="right">SourceSeq</TrackingResultTable>
                                            <TrackingResultTable align="right">Identifier</TrackingResultTable>
                                            <TrackingResultTable align="right">Source country</TrackingResultTable>
                                            <TrackingResultTable align="right">Isolated year</TrackingResultTable>
                                            <TrackingResultTable align="right">Serogroup type</TrackingResultTable>
                                            <TrackingResultTable align="right">ST</TrackingResultTable>
                                            <TrackingResultTable align="right">Clonal complex</TrackingResultTable>
                                            <TrackingResultTable align="right">PorA</TrackingResultTable>
                                            <TrackingResultTable align="right">PorB</TrackingResultTable>
                                            <TrackingResultTable align="right">fHbp</TrackingResultTable>
                                            <TrackingResultTable align="right">NHBA</TrackingResultTable>
                                            <TrackingResultTable align="right">NadA</TrackingResultTable>
                                            <TrackingResultTable align="right">FetA</TrackingResultTable>
                                            <TrackingResultTable align="right">BAST</TrackingResultTable>
                                            <TrackingResultTable align="right">ResistanceGene</TrackingResultTable>
                                            <TrackingResultTable align="right">No. void loci</TrackingResultTable>
										</TableRow>
									</TableHead>
									<TableBody>
										{trackResult.map(row => (
											<TableRow className={classes.row} key={row.BioSample}>
												<TrackingResultTable align="right">{row["Distance (loci)"]}</TrackingResultTable>
												<TrackingResultTable align="right">
													<a href={'https://www.ncbi.nlm.nih.gov/biosample/?term='+row.BioSample}
													target="_blank">
														{row.BioSample}
													</a>
												</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceSeq}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.Identifier}</TrackingResultTable>
												<TrackingResultTable align="right">{row.SourceCountry}</TrackingResultTable>
												<TrackingResultTable align="right">{row.IsolatYear}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Serogroup_type}</TrackingResultTable>
												<TrackingResultTable align="right">{row.ST}</TrackingResultTable>
												<TrackingResultTable align="right">{row["Clonal complex"]}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.PorA}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.PorB}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.fHbp}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.NHBA}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.NadA}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.FetA}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.BAST}</TrackingResultTable>
                                                <TrackingResultTable align="right">{row.ResistanceGene}</TrackingResultTable>
												<TrackingResultTable align="right">{row.Void_loci}</TrackingResultTable>
											</TableRow>
										))}
									</TableBody>
								</Table>
							</Scrollbars>
							<br />
						</Paper>
					</div>
					<br />
					<div className={classes.displayCenter}>
						<Link to="/cgMLST/tracking" style={{ textDecoration:'none' }}>
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

export default withStyles(styles)(TrackingResult);