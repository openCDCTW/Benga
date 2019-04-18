import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
//Option's outlook
import PropTypes from 'prop-types';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
//
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import blue from '@material-ui/core/colors/blue';
//search
import TextField from '@material-ui/core/TextField';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import SearchIcon from '@material-ui/icons/Search';
import Typography from '@material-ui/core/Typography';
import OutlinedInput from '@material-ui/core/OutlinedInput';

const styles = theme => ({
    buttoncss:{
        color: theme.palette.getContrastText(blue[600]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
    },
    selectcss: {
	    display: 'flex',
	    flexWrap: 'wrap',
	},
	formControl: {
	    margin: theme.spacing.unit,
	    minWidth: 240,
	},
	selectEmpty: {
	    marginTop: theme.spacing.unit * 2,
	},
        textField:{
        marginLeft: '19px',
        marginTop: '25px',
        width:'22.5%'
    },
    container:{
        display:'flex',
        flexWrap:'wrap',
    },
    numberField:{
        marginLeft: '19px',
        marginTop: '20px',
        marginRight: '10px',
        width:'110px'
    },
    countrySelect:{
        marginLeft: '19px',
        marginTop: '25px',
        width: '22.5%',
    }
    
})


class Tracking extends React.Component {

	constructor(props) {

		super(props);

        let nowYear = new Date();

		this.state = {
            allele_db:"Vibrio_cholerae",
            profile_db:"Vibrio_cholerae",
            yearError: false,
            country: '',
            labelWidth: 0,
            nowYear: nowYear.getFullYear(),
		};

		this.djsConfig = {
            dictDefaultMessage:"Drag a cgMLST profile here",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 1,
            init:function(){
                this.on("addedfile", function(on_load_header_data){
                    
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("profile_db", "Vibrio_cholerae");
                });
                this.on("success", function(file,response){
                    file._removeLink.remove();
                    delete file._removeLink;
                    window.trackingID = response.id;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/tracking/profile/'
        };

        this.dropzone = null;
	}

    componentDidMount() {
        // this.setState({
        //   labelWidth: ReactDOM.findDOMNode(this.InputLabelRef).offsetWidth,
        // });
    }

    handlePost() {
    
	    let fileCheck = this.dropzone.files.length;

	    if(fileCheck != 1){
	        alert('Please upload only 1 file');
	        return ;
	    }
	    this.dropzone.processQueue();
        this.props.history.push("/tracking_result")
	}

	remove(){
	    this.dropzone.removeAllFiles();
	}

	select_handleChange(event){
        if( event.target.value == 'Vibrio_cholerae'){
            this.setState(state => ({ 
                [event.target.name]: event.target.value,
                profile_db:"vibrio-profiles",
            }));
        }
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

    serotypeHandleChange(event){
        this.setState(state => ({ serotype: this.serotype.value }));
    }

    search(){
        console.log(this.state.bioSampleID);
        console.log(this.state.strain);
        console.log(Number.parseInt(this.state.yearFrom, 10));
        console.log(Number.parseInt(this.state.yearTo, 10));
        console.log(this.state.serotype);

        if( this.state.yearFrom > this.state.yearTo ){
            this.setState(state => ({ yearError: true }));
            alert('Input year error');
        }
    }

	render() {

		const config = this.componentConfig;
		const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        return (
            <div>
                <br />
                <br />
                <div>
                    <div style={{ float:'left', marginLeft:'10px' }}>
                        <form className={classes.selectcss} autoComplete="off">
                            <FormControl required className={classes.formControl} >
                              <InputLabel htmlFor="database-required">Database</InputLabel>
                                <Select
                                  value={this.state.allele_db}
                                  onChange={this.select_handleChange.bind(this)}
                                  name="allele_db"
                                  inputProps={{
                                    id: 'database-required',
                                  }}
                                  className={classes.selectEmpty}
                                  >
                                  <MenuItem value={'Vibrio_cholerae'}>Vibrio cholerae</MenuItem>
                                </Select>
                              <FormHelperText>Required</FormHelperText>
                            </FormControl>
                        </form>
                    </div>
                    <div style={{ float:'right', marginTop:'35px', marginRight:'25px' }}>
                        <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                                Remove all files
                                &nbsp;&nbsp;
                                <DeleteIcon />
                        </Button>
                    </div>
                </div>
                <br />
                <br />
                <br />
                <br />
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                </div>
                <br />
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" className ={classes.buttoncss} 
                     onClick={this.handlePost.bind(this)}>
                        Submit
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                </div>
                <br />
                <br />
                <br />
                <br />
            </div>
        );
    }
}

export default withStyles(styles)(Tracking);



                // <div>
                //     <Paper>
                //         <div>
                //             <br />
                //             <Typography style={{ marginLeft:'20px', fontSize:'25px', fontWeight:'500' }}>
                //                 Search database
                //             </Typography>
                //             <form className={classes.container}>
                //                 <TextField
                //                     inputRef={ID => this.biosample = ID}
                //                     label="BioSample ID"
                //                     type="search"
                //                     placeholder="Input BioSample ID here"
                //                     onChange={this.bioSampleHandleChange.bind(this)}
                //                     className={classes.textField}
                //                     onKeyPress={this._onKeyPress}
                //                     margin="normal"
                //                     variant="outlined"
                //                 />
                //                 <TextField
                //                     inputRef={ID => this.strain = ID}
                //                     label="Strain ID"
                //                     type="search"
                //                     placeholder="Input strain ID here"
                //                     onChange={this.strainHandleChange.bind(this)}
                //                     className={classes.textField}
                //                     onKeyPress={this._onKeyPress}
                //                     margin="normal"
                //                     variant="outlined"
                //                 />
                //                 <TextField
                //                     inputRef={ID => this.serotype = ID}
                //                     label="ST"
                //                     type="search"
                //                     placeholder="Input serotype here"
                //                     onChange={this.serotypeHandleChange.bind(this)}
                //                     className={classes.textField}
                //                     onKeyPress={this._onKeyPress}
                //                     margin="normal"
                //                     variant="outlined"
                //                 />
                //                 <FormControl variant="outlined" className={classes.countrySelect}>
                //                     <InputLabel
                //                         ref={ref => {
                //                           this.InputLabelRef = ref;
                //                         }}
                //                     >
                //                     Country
                //                     </InputLabel>
                //                     <Select
                //                     value={this.state.country}
                //                     onChange={this.countryHandleChange.bind(this)}
                //                     name="country"
                //                     input={
                //                             <OutlinedInput
                //                             labelWidth={this.state.labelWidth}
                //                             name="country"
                //                             />
                //                         }
                //                     MenuProps={{
                //                         PaperProps: {
                //                             style: {
                //                                 maxHeight: 49 * 4.5 + 8,
                //                             },
                //                           },
                //                       }}
                //                     >
                //                         <MenuItem value=''>Any</MenuItem>
                //                         <MenuItem value='Taiwan'>Taiwan</MenuItem>
                //                         <MenuItem value='China'>China</MenuItem>
                //                         <MenuItem value='Russia'>Russia</MenuItem>
                //                         <MenuItem value='Ukraine'>Ukraine</MenuItem>
                //                         <MenuItem value='USA'>USA</MenuItem>
                //                     </Select>
                //                 </FormControl>
                //             </form>
                //         </div>
                //         <div>
                //             <form className={classes.container}>
                //                 <TextField
                //                     label="Year"
                //                     inputRef={year => this.yearFrom = year}
                //                     type="number"
                //                     onChange={this.yearFromHandleChange.bind(this)}
                //                     className={classes.numberField}
                //                     inputProps={{ 
                //                         style:{ textAlign: 'center'},
                //                     }}
                //                     margin="normal"
                //                     variant="outlined"
                //                     error={this.state.yearError}
                //                 />
                //                 &nbsp;&nbsp;
                //                 <Typography style={{ fontSize:'18px', display:'flex', justifyContent:'center', 
                //                 alignItems:'center'}}>to</Typography>
                //                 &nbsp;&nbsp;
                //                 <TextField
                //                     label="Year"
                //                     inputRef={year => this.yearTo = year}
                //                     type="number"
                //                     onChange={this.yearToHandleChange.bind(this)}
                //                     className={classes.numberField}
                //                     inputProps={{ 
                //                         style:{ textAlign: 'center'},
                //                     }}
                //                     margin="normal"
                //                     variant="outlined"
                //                     error={this.state.yearError}
                //                 />
                //             </form>
                //         </div>
                //         <br />
                //         <div style={{ width:'97%', display:'flex', justifyContent:'flex-end', 
                //                 alignItems:'flex-end'}}>
                //             <Button variant="contained" color="default" onClick={this.search.bind(this)}>
                //                 Search
                //                 &nbsp;&nbsp;
                //                 <SearchIcon />
                //             </Button>
                //         </div>
                //         <br />
                //     </Paper>
                //     <br />
                //     <br />
                //     <br />
                //     <br />
                // </div>