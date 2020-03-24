import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route, Redirect } from "react-router-dom";
import { withStyles } from '@material-ui/core/styles';
import Drawer from '@material-ui/core/Drawer';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import Divider from '@material-ui/core/Divider';
import FunctionLists from './FunctionLists.jsx';

import About from './Functions/About.jsx';
import Profiling from './Functions/Profiling.jsx';
import ProfilingResult from './Functions/ProfilingResult.jsx';
import Tracking from './Functions/Tracking.jsx';
import TrackingResult from './Functions/TrackingResult.jsx';
import Clustering from './Functions/Clustering.jsx';
import ClusteringResult from './Functions/ClusteringResult.jsx';
import TrackingSearch from './Functions/TrackingSearch.jsx';
import QuerybyID from './Functions/QuerybyID.jsx';
import Footer from './Footer.jsx';

const drawerWidth = 300;
const styles = theme => ({
    root: {
        margin: '20px',
        display: 'flex',
        flexWrap: 'wrap',
    },
    drawer: {
        width: drawerWidth,
        flexShrink: 0,
    },
    drawerPaper: {
        width: drawerWidth,
    },
    formControl: {
        marginLeft: theme.spacing(1),
        minWidth: 240,
    },
    selectEmpty: {
        marginTop: theme.spacing(2),
    },
    content: {
        marginLeft: '300px',
        height: '100vh',
        width: '80vw',
    },
    toolbar: theme.mixins.toolbar,
});

class PageContent extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            database: 'Vibrio_cholerae',
            toProfiling: false,
        };
        window.database = 'Vibrio_cholerae';
        window.occurrence = 95;
    }

    handleChange(event){
        this.setState(state => ({ [event.target.name]: event.target.value,
                                   toProfiling: true,
                                }));
        window.database = event.target.value;
        if(event.target.value == 'Salmonella_enterica'){
            window.occurrence = 92;
        }else{
            window.occurrence = 95;
        }
    };

    render() {
        const { classes } = this.props;

        return(
            <Router>
                <div>
                    <Drawer
                    className={classes.drawer}
                    variant="permanent"
                    classes={{
                        paper: classes.drawerPaper,
                    }}
                    >
                        <div className={classes.toolbar} />
                        <form className={classes.root} autoComplete="off">
                            <FormControl required className={classes.formControl}>
                                <InputLabel htmlFor="database-required">Species</InputLabel>
                                <Select
                                    value={this.state.database}
                                    onChange={this.handleChange.bind(this)}
                                    name="database"
                                    inputProps={{
                                        id: 'database-required',
                                    }}
                                    className={classes.selectEmpty}>
                                        <MenuItem value={'Vibrio_cholerae'}>Vibrio cholerae</MenuItem>
                                        <MenuItem value={'Neisseria_meningitidis'}>Neisseria meningitidis</MenuItem>
                                        <MenuItem value={'Salmonella_enterica'}>Salmonella enterica</MenuItem>
                                </Select>
                                <FormHelperText>Required</FormHelperText>
                            </FormControl>
                        </form>
                        <Divider />
                        <FunctionLists species={this.state.database} />
                    </Drawer>
                    <main className={classes.content}>
                        <div className={classes.toolbar} />
                        { this.state.toProfiling == true ? <Redirect to='/cgMLST/' /> : null }
                        <Route path="/cgMLST/about" exact component={() => <About />} />
                        <Route path="/cgMLST/" exact component={() => <Profiling />} />
                        <Route path="/cgMLST/profiling" exact component={() => <Profiling />} />
                        <Route path="/cgMLST/profilingResult" exact component={() => <ProfilingResult />} />
                        <Route path="/cgMLST/tracking" exact component={() => <Tracking />} />
                        <Route path="/cgMLST/trackingResult" exact component={() => <TrackingResult />} />
                        <Route path="/cgMLST/Search" exact component={() => <TrackingSearch />} />
                        <Route path="/cgMLST/clustering" component={() => <Clustering />} />
                        <Route path="/cgMLST/clusteringResult" component={() => <ClusteringResult />} />
                        <Route path="/cgMLST/queryByID" component={() => <QuerybyID />} />
                        <Footer />
                    </main>
                </div>
            </Router>
        );
    }
}

export default withStyles(styles)(PageContent);