import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TfDataComponent } from './tf-data.component';

describe('TfDataComponent', () => {
  let component: TfDataComponent;
  let fixture: ComponentFixture<TfDataComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TfDataComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(TfDataComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
